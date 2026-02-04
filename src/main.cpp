#include <netcdf.h>
#include <cstdio>
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include <string>
#include <sstream>

#include "imgui.h"
#include "implot.h"

#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"

#include <GLFW/glfw3.h>

// If GLFW fails, this prints a readable error
static void glfw_error_callback(int error, const char* description) {
    std::fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

// helper: split csv line by comma
static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) out.push_back(cell);
    return out;
}

static bool starts_with(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

struct NcHandle {
    int ncid = -1;
    int var_temp = -1;
    int var_lat = -1;
    int var_lon = -1;
    int var_time = -1;

    int dim_time = -1;
    int dim_lat = -1;
    int dim_lon = -1;

    size_t ntime = 0, nlat = 0, nlon = 0;
    bool ok = false;
};

static bool nc_init(NcHandle& h, const char* path) {
    if (nc_open(path, NC_NOWRITE, &h.ncid) != NC_NOERR) return false;

    if (nc_inq_varid(h.ncid, "tempanomaly", &h.var_temp) != NC_NOERR) return false;
    if (nc_inq_varid(h.ncid, "lat", &h.var_lat) != NC_NOERR) return false;
    if (nc_inq_varid(h.ncid, "lon", &h.var_lon) != NC_NOERR) return false;
    if (nc_inq_varid(h.ncid, "time", &h.var_time) != NC_NOERR) return false;

    // dims from variables (more robust than hardcoding dim names)
    int ndims = 0;
    int dimids[NC_MAX_VAR_DIMS];

    if (nc_inq_var(h.ncid, h.var_temp, nullptr, nullptr, &ndims, dimids, nullptr) != NC_NOERR) return false;
    // expected order: time, lat, lon
    h.dim_time = dimids[0];
    h.dim_lat  = dimids[1];
    h.dim_lon  = dimids[2];

    if (nc_inq_dimlen(h.ncid, h.dim_time, &h.ntime) != NC_NOERR) return false;
    if (nc_inq_dimlen(h.ncid, h.dim_lat,  &h.nlat)  != NC_NOERR) return false;
    if (nc_inq_dimlen(h.ncid, h.dim_lon,  &h.nlon)  != NC_NOERR) return false;

    h.ok = true;
    return true;
}

static void nc_shutdown(NcHandle& h) {
    if (h.ncid != -1) nc_close(h.ncid);
    h = NcHandle{};
}


static bool nc_read_slab(const NcHandle& h, int time_index, std::vector<short>& out) {
    if (!h.ok) return false;
    out.resize(h.nlat * h.nlon);

    size_t start[3] = { (size_t)time_index, 0, 0 };
    size_t count[3] = { 1, h.nlat, h.nlon };

    int err = nc_get_vara_short(h.ncid, h.var_temp, start, count, out.data());
    return err == NC_NOERR;
}


static bool nc_read_cell_series(const NcHandle& h, int ilat, int ilon, std::vector<short>& out) {
    if (!h.ok) return false;
    out.resize(h.ntime);

    size_t start[3] = { 0, (size_t)ilat, (size_t)ilon };
    size_t count[3] = { h.ntime, 1, 1 };

    int err = nc_get_vara_short(h.ncid, h.var_temp, start, count, out.data());
    return err == NC_NOERR;
}


static void build_lat_weights(const std::vector<float>& lat, std::vector<double>& w_row) {
    w_row.resize(lat.size());
    for (size_t i = 0; i < lat.size(); ++i) {
        double rad = lat[i] * 3.141592653589793 / 180.0;
        w_row[i] = std::cos(rad);
        if (w_row[i] < 0) w_row[i] = 0; // just in case of weirdness
    }
}



struct ClimateSeries {
    std::vector<double> years;
    std::vector<double> jd_anomaly;
};

static ClimateSeries load_gistemp_global_jd(const std::string& path){
    ClimateSeries series;

    std::ifstream file(path);
    if (!file.is_open()) {
        std::fprintf(stderr, "Failed to open data file: %s\n", path.c_str());
        return series;
    }

    std::string line;
    std::vector<std::string> header;
    int jd_index = -1;

    while (std::getline(file, line)) {
        if (starts_with(line, "Year,")) {
            header = split_csv_line(line);
            for (int i = 0; i < (int)header.size(); i++) {
                if (header[i] == "J-D") {
                    jd_index = i;
                    break;
                }
            }
            break;
        }
    }

    if (jd_index == -1) {
        std::fprintf(stderr, "J-D column not found in data file: %s\n", path.c_str());
        return series;
    }

    // Read data lines
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        auto cells = split_csv_line(line);
        if (cells.size() <= (size_t)jd_index) continue;

        if (cells[0].empty() || !(cells[0][0] >= '0' && cells[0][0] <= '9')) break;

        // Parse year
        int year = 0;
        try {
            year = std::stoi(cells[0]);
        } catch (...) {
            continue;
        }

        // Parse anomaly (skip missing values like "***")
        const std::string& jd = cells[jd_index];
        if (jd.empty() || jd.find('*') != std::string::npos) continue;

        try{
            double val = std::stod(jd);
            series.years.push_back((double)year);
            series.jd_anomaly.push_back(val);
        }
        catch (...) {
            continue;
        }
    }
    return series;

}

static int find_nearest_year_index(const std::vector<double>& years, double year_query) {
    if (years.empty()) return -1;

    int best_i = 0;
    double best_d = std::abs(years[0] - year_query);

    for (int i = 1; i < (int)years.size(); i++) {
        double d = std::abs(years[i] - year_query);
        if (d < best_d) {
            best_d = d;
            best_i = i;
        }
    }
    return best_i;
}

struct Grid2D {
    int nlat = 0;
    int nlon = 0;
    std::vector<float> lat;
    std::vector<float> lon;
    std::vector<float> values; // size = nlat * nlon
};


static bool nc_ok(int status, const char* msg){
    if (status != NC_NOERR) {
        std::fprintf(stderr, "netCDF error (%s): %s\n", msg, nc_strerror(status));
        return false;
    }
    return true;
}


static bool load_gistemp_grid_t_from_handle(const NcHandle& nc, int time_index, Grid2D& out) {
    if (!nc.ok) return false;
    if (time_index < 0 || (size_t)time_index >= nc.ntime) return false;

    out.nlat = (int)nc.nlat;
    out.nlon = (int)nc.nlon;

    out.lat.resize(nc.nlat);
    out.lon.resize(nc.nlon);

    if (nc_get_var_float(nc.ncid, nc.var_lat, out.lat.data()) != NC_NOERR) return false;
    if (nc_get_var_float(nc.ncid, nc.var_lon, out.lon.data()) != NC_NOERR) return false;

    std::vector<short> raw;
    if (!nc_read_slab(nc, time_index, raw)) return false;

    const short fill = 32767;
    const float scale = 0.01f;

    out.values.resize(nc.nlat * nc.nlon);
    for (size_t i = 0; i < raw.size(); ++i) {
        out.values[i] = (raw[i] == fill) ? std::numeric_limits<float>::quiet_NaN()
                                         : (float)raw[i] * scale;
    }
    return true;
}


static unsigned char clamp_u8(int x){
    if (x < 0) return 0;
    if (x > 255) return 255;
    return (unsigned char)x;
}

struct TimeAxis{
    std::vector<int> days_since_1800; //size ntime
};


static bool load_time_axis_from_handle(const NcHandle& nc, TimeAxis& out) {
    if (!nc.ok) return false;

    out.days_since_1800.resize(nc.ntime);
    int err = nc_get_var_int(nc.ncid, nc.var_time, out.days_since_1800.data());
    return err == NC_NOERR;
}



static bool is_leap(int y) {
    if (y % 400 == 0) return true;
    if (y % 100 == 0) return false;
    return (y % 4 == 0);
}


static void days_since_1800_to_ym(int days, int& year, int& month) {
    // Start at 1800-01-01
    year = 1800;
    month = 1;

    int d = days;
    while (true) {
        int days_in_year = is_leap(year) ? 366 : 365;
        if (d >= days_in_year) { d -= days_in_year; year++; }
        else break;
    }

    static const int mdays_norm[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
    for (int m = 1; m <= 12; m++) {
        int dim = mdays_norm[m-1];
        if (m == 2 && is_leap(year)) dim = 29;
        if (d >= dim) { d -= dim; month++; }
        else break;
    }
    // month is now 1..12
}


static const char* month_name(int m) {
    static const char* names[12] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
    if (m < 1 || m > 12) return "???";
    return names[m-1];
}



static void anomaly_to_rgb(float v, float vmin, float vmax,
                           unsigned char& r, unsigned char& g, unsigned char& b) {
    if (std::isnan(v)) { r = g = b = 0; return; } // missing -> black

    if (vmax <= vmin) vmax = vmin + 1e-6f;

    float t = (v - vmin) / (vmax - vmin); // vmin..vmax -> 0..1
    if (t < 0) t = 0;
    if (t > 1) t = 1;

    // blue -> white -> red
    if (t < 0.5f) {
        float u = t / 0.5f;
        r = clamp_u8((int)(255.0f * u));
        g = clamp_u8((int)(255.0f * u));
        b = 255;
    } else {
        float u = (t - 0.5f) / 0.5f;
        r = 255;
        g = clamp_u8((int)(255.0f * (1.0f - u)));
        b = clamp_u8((int)(255.0f * (1.0f - u)));
    }
}




static void draw_anomaly_legend(float vmin=-2.0f, float vmax=2.0f, ImVec2 size=ImVec2(260, 16)) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    ImVec2 p0 = ImGui::GetCursorScreenPos();
    ImVec2 p1(p0.x + size.x, p0.y + size.y);

    const int steps = 200;
    for (int i = 0; i < steps; i++) {
        float t0 = (float)i / (float)steps;
        float t1 = (float)(i + 1) / (float)steps;

        float v0 = vmin + (vmax - vmin) * t0;

        unsigned char r,g,b;
        anomaly_to_rgb(v0, vmin, vmax, r, g, b);

        ImU32 col = IM_COL32(r, g, b, 255);

        float x0 = p0.x + size.x * t0;
        float x1 = p0.x + size.x * t1;

        dl->AddRectFilled(ImVec2(x0, p0.y), ImVec2(x1, p1.y), col);
    }

    dl->AddRect(p0, p1, IM_COL32(255,255,255,120));

    ImGui::Dummy(size);
    ImGui::Text("Scale: %.1fC .. %.1fC", vmin, vmax);
}


static GLuint create_texture_rgba(int w, int h, const std::vector<unsigned char>& rgba) {
    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgba.data());
    return tex;
}

static void update_texture_rgba(GLuint tex, int w, int h, const std::vector<unsigned char>& rgba) {
    glBindTexture(GL_TEXTURE_2D, tex);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, rgba.data());
}

struct CellSeries{
    int ilat = -1;
    int ilon = -1;
    std::vector<double> x_year;
    std::vector<double> y_anom;
    bool valid = false;
};


static bool load_cell_timeseries_from_handle(
    const NcHandle& nc,
    int ilat,
    int ilon,
    const TimeAxis& t,
    CellSeries& out
) {
    if (!nc.ok) return false;

    if (ilat < 0 || ilon < 0 || (size_t)ilat >= nc.nlat || (size_t)ilon >= nc.nlon) return false;
    if (t.days_since_1800.size() != nc.ntime) return false;

    std::vector<short> raw;
    if (!nc_read_cell_series(nc, ilat, ilon, raw)) return false;

    const short fill = 32767;
    const double scale = 0.01;

    out.ilat = ilat;
    out.ilon = ilon;
    out.x_year.resize(nc.ntime);
    out.y_anom.resize(nc.ntime);

    for (size_t i = 0; i < nc.ntime; i++) {
        int y = 0, m = 0;
        days_since_1800_to_ym(t.days_since_1800[i], y, m);
        out.x_year[i] = (double)y + ((double)m - 0.5) / 12.0;

        out.y_anom[i] = (raw[i] == fill) ? std::numeric_limits<double>::quiet_NaN()
                                         : (double)raw[i] * scale;
    }

    out.valid = true;
    return true;
}


int main() {

    ClimateSeries climate_data = load_gistemp_global_jd("data/global_temp.csv");

    if (climate_data.years.empty()) {
        std::fprintf(stderr, "No climate data loaded.\n");
        return 1;
    }
    // 1) Initialize GLFW (windowing + input)
    glfwSetErrorCallback(glfw_error_callback);

    if (!glfwInit()) {
        std::fprintf(stderr, "Failed to init GLFW\n");
        return 1;
    }

    // 2) Create an OpenGL context (OpenGL 2.1 for simplicity)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

    GLFWwindow* window = glfwCreateWindow(1000, 700, "Climate Lab", nullptr, nullptr);
    if (!window) {
        std::fprintf(stderr, "Failed to create window\n");
        glfwTerminate();
        return 1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // VSync: caps FPS to your display refresh rate

    // 3) Create ImGui + ImPlot contexts (their global state)
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();

    // 4) Choose a default UI style
    ImGui::StyleColorsDark();

    // 5) Connect ImGui to GLFW + OpenGL2 (“backend” glue)
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();

    Grid2D grid;
    const char* nc_path = "data/gistemp1200_GHCNv4_ERSSTv5.nc";
    NcHandle nc;
    if (!nc_init(nc, nc_path)) {
        std::fprintf(stderr, "Failed to init NetCDF handle.\n");
        return 1;
    }
    std::fprintf(stderr, "NetCDF dims: time=%zu lat=%zu lon=%zu\n", nc.ntime, nc.nlat, nc.nlon);

    TimeAxis t;
    if (!load_time_axis_from_handle(nc, t)) {
        std::fprintf(stderr, "Failed to load time axis.\n");
        return 1;
    }

    int time_index = 0;
    int year = 0, month = 0;
    days_since_1800_to_ym(t.days_since_1800[time_index], year, month);

    // Pick a starting time index 
    if (!load_gistemp_grid_t_from_handle(nc, time_index, grid)) {
        std::fprintf(stderr, "Failed to load grid from netCDF.\n");
        return 1;
    }


    // Build RGBA image
    std::vector<unsigned char> rgba(grid.nlat * grid.nlon * 4);

    float vmin = -2.0f;
    float vmax = 2.0f;
    
    auto rebuild_rgba = [&]() {
        for (int ilat = 0; ilat < grid.nlat; ilat++) {
            for (int ilon = 0; ilon < grid.nlon; ilon++) {
                int src_ilat = ilat;
                int dst_ilat = (grid.nlat - 1 - ilat);

                float v = grid.values[src_ilat * grid.nlon + ilon];

                unsigned char r, g, b;
                anomaly_to_rgb(v, vmin, vmax, r, g, b);

                int idx = (dst_ilat * grid.nlon + ilon) * 4;
                rgba[idx + 0] = r;
                rgba[idx + 1] = g;
                rgba[idx + 2] = b;
                rgba[idx + 3] = 255;
            }
        }
    };

    rebuild_rgba();
    GLuint tex = 0;
    tex = create_texture_rgba(grid.nlon, grid.nlat, rgba);

    

    // 7) Main loop: handle input → build UI → render → repeat
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // Start a new UI frame
        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Build the UI (widgets + plots)
        ImGui::Begin("Dashboard");
        ImGui::Text("Climate Lab v0.1 (UI online)");
        ImGui::Separator();

        if (ImPlot::BeginPlot("Global Temperature Anomaly (NASA GISTEMP)")) {
            ImPlot::PlotLine("J-D anomaly", climate_data.years.data(), climate_data.jd_anomaly.data(),
                 (int)climate_data.years.size());
            
            //hover interaction
            if (ImPlot::IsPlotHovered()) {
                ImPlotPoint mouse = ImPlot::GetPlotMousePos(); // in plot units (Year, C)

                int idx = find_nearest_year_index(climate_data.years, mouse.x);
                if (idx >= 0) {
                    double year = climate_data.years[idx];
                    double anom = climate_data.jd_anomaly[idx];

                    ImPlotPoint p_top(year, ImPlot::GetPlotLimits().Y.Max);
                    ImPlotPoint p_bot(year, ImPlot::GetPlotLimits().Y.Min);

                    ImVec2 pix_top = ImPlot::PlotToPixels(p_top);
                    ImVec2 pix_bot = ImPlot::PlotToPixels(p_bot);

                    // Draw
                    ImDrawList* draw_list = ImGui::GetWindowDrawList();
                    draw_list->AddLine(pix_top, pix_bot, IM_COL32(255, 255, 255, 255), 1.0f);

                    // Highlight point
                    ImPlot::PlotScatter("##hover_point", &year, &anom, 1);

                    ImGui::BeginTooltip();
                    ImGui::Text("Year: %.0f", year);
                    ImGui::Text("J-D anomaly: %.3f C", anom);
                    ImGui::EndTooltip();
                }
            }
            ImPlot::EndPlot();

        }

        ImGui::End();

        ImGui::Begin("Map View");


        static int sel_ilat = -1;
        static int sel_ilon = -1;
        static CellSeries cell_ts;
        static bool focus_cell_window = false;

        ImGui::Text("Time: %s %d (index %d)", month_name(month), year, time_index);
        ImGui::SliderInt("Time##time_index", &time_index, 0, (int)t.days_since_1800.size() -1);

        static int last_time_index = -1;
        if (time_index != last_time_index) {
            last_time_index = time_index;
            

            days_since_1800_to_ym(t.days_since_1800[time_index], year, month);
            if (load_gistemp_grid_t_from_handle(nc, time_index, grid)) {
                rebuild_rgba();
                update_texture_rgba(tex, grid.nlon, grid.nlat, rgba);
            }
        }

        float avail_w = ImGui::GetContentRegionAvail().x;
        float aspect = (float)grid.nlat / (float)grid.nlon;
        ImVec2 size(avail_w, avail_w * aspect);
        ImGui::Text("GISTEMP gridded anomaly (time index %d)", time_index);
        ImGui::Image((ImTextureID)(intptr_t)tex, size);



        
        ImVec2 img_min = ImGui::GetItemRectMin();
        ImVec2 img_max = ImGui::GetItemRectMax();
        float img_w = img_max.x - img_min.x;
        float img_h = img_max.y - img_min.y;

        float cell_w = img_w / (float)grid.nlon;
        float cell_h = img_h / (float)grid.nlat;

        if (sel_ilat >= 0 && sel_ilon >= 0){
            int py = (grid.nlat -1 - sel_ilat);
            int px = sel_ilon;
            ImVec2 p0(img_min.x + px * cell_w, img_min.y + py * cell_h);
            ImVec2 p1(p0.x + cell_w, p0.y + cell_h);
            ImGui::GetWindowDrawList()->AddRect(p0, p1, IM_COL32(0,255,0,255), 0.0f, 0, 2.0f);
        }
        
        

        if (ImGui::IsItemHovered()) {

            // Mouse position in screen coords
            ImVec2 mouse = ImGui::GetIO().MousePos;

            // Normalized mouse position inside the image (0..1)
            float u = (mouse.x - img_min.x) / (img_max.x - img_min.x);
            float v = (mouse.y - img_min.y) / (img_max.y - img_min.y);

            // Clamp to [0,1] in case of tiny overshoot
            if (u < 0) u = 0; if (u > 1) u = 1;
            if (v < 0) v = 0; if (v > 1) v = 1;

            // Convert to pixel indices in the *rendered image*
            int px = (int)std::floor(u * grid.nlon);  // 0..nlon-1
            int py = (int)std::floor(v * grid.nlat);  // 0..nlat-1

            if (px >= grid.nlon) px = grid.nlon - 1;
            if (py >= grid.nlat) py = grid.nlat - 1;

            // Because we flipped latitude when building the texture,
            // the rendered image row py corresponds to source latitude index:
            int ilon = px;
            int ilat_src = (grid.nlat - 1 - py);

            float lat = grid.lat[ilat_src];
            float lon = grid.lon[ilon];
            float val = grid.values[ilat_src * grid.nlon + ilon];
            float lon_disp = lon;
            if (lon_disp > 180.0f) lon_disp -= 360.0f;


            if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && !ImGui::IsAnyItemActive()) {
                sel_ilat = ilat_src;
                sel_ilon = ilon;

                if (!load_cell_timeseries_from_handle(nc, sel_ilat, sel_ilon, t, cell_ts)) {
                    std::fprintf(stderr, "Failed to load cell time series.\n");
                    cell_ts.valid = false;
                }

                focus_cell_window = true;
        }


        ImGui::BeginTooltip();
        ImGui::Text("Time: %s %d", month_name(month), year);
        ImGui::Separator();
        ImGui::Text("lat: %.2f", lat);
        ImGui::Text("lon: %.2f", lon_disp);

        if (std::isnan(val)) {
            ImGui::Text("anomaly: (missing)");
        } else {
            ImGui::Text("anomaly: %.2f C", val);
        }

        ImGui::Text("cell: ilat=%d ilon=%d", ilat_src, ilon);
        ImGui::EndTooltip();
        }


        if (ImGui::Button("Reset scale")) {
        vmin = -2.0f;
        vmax =  2.0f;
        rebuild_rgba();
        update_texture_rgba(tex, grid.nlon, grid.nlat, rgba);
        }
        
        bool scale_changed = false;
        scale_changed |= ImGui::SliderFloat("Vmin", &vmin, -5.0f, 0.0f, "%.2f C");
        scale_changed |= ImGui::SliderFloat("Vmax", &vmax,  0.0f, 5.0f, "%.2f C");

        if (vmax <= vmin) vmax = vmin + 0.01f;

        // legend full width
        float legend_w = ImGui::GetContentRegionAvail().x;
        draw_anomaly_legend(vmin, vmax, ImVec2(legend_w, 16));

        if (scale_changed) {
            rebuild_rgba();
            update_texture_rgba(tex, grid.nlon, grid.nlat, rgba);
        }
        

        ImGui::End();

        if (focus_cell_window) {
            ImGui::SetNextWindowFocus();
            focus_cell_window = false;
        }

        ImGui::Begin("Cell Time Series");

        if (sel_ilat >= 0 && sel_ilon >= 0) {
            float lat_sel = grid.lat[sel_ilat];
            float lon_sel = grid.lon[sel_ilon];
            float lon_disp = lon_sel;
            if (lon_disp > 180.0f) lon_disp -= 360.0f;

            ImGui::Text("Selected cell: ilat=%d ilon=%d | lat=%.2f lon=%.2f", sel_ilat, sel_ilon, lat_sel, lon_disp);
            float cur = grid.values[sel_ilat * grid.nlon + sel_ilon];
            if (!std::isnan(cur))
                ImGui::Text("Current month: %.2f °C", cur);
            else
                ImGui::Text("Current month: (missing)");


            if (cell_ts.valid && !cell_ts.x_year.empty()){
                if (ImPlot::BeginPlot("Temp anomaly at selected cell")) {
                    ImPlot::PlotLine("anomaly", cell_ts.x_year.data(), cell_ts.y_anom.data(), (int)cell_ts.x_year.size());
                    ImPlot::EndPlot();
                }
            } else {
                ImGui::Text("No valid data for this cell.");
            }
        } else {
            ImGui::Text("No cell selected. Click on the map to select a cell.");
        }




        ImGui::Text("Index %d / %d", time_index, (int)t.days_since_1800.size() - 1);

        ImGui::Text("Monthly mean anomaly (Base: 1951 - 1980)");
        ImGui::Text("lat[0]=%.2f lat[last]=%.2f", grid.lat.front(), grid.lat.back());
        ImGui::Text("lon[0]=%.2f lon[last]=%.2f", grid.lon.front(), grid.lon.back());


        
        ImGui::End();

        // Render UI to OpenGL
        ImGui::Render();

        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // 8) Cleanup: shutdown backends + destroy contexts
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    nc_shutdown(nc);
    return 0;
}
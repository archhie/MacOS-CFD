#include "colormaps.hpp"
#include <algorithm>
#include <cmath>

namespace Colormaps {

void viridis(float t, float& r, float& g, float& b) {
    // Viridis colormap implementation
    const float c0[3] = {0.2777273272234177f, 0.005407344544966578f, 0.3340998053353061f};
    const float c1[3] = {0.1050930431085774f, 1.404613529898575f, 1.384590162594685f};
    const float c2[3] = {-0.3308618287255563f, 0.214847559468213f, 0.09509516302823659f};
    const float c3[3] = {-4.634230498983486f, -5.799100973351585f, -19.33244095627987f};
    const float c4[3] = {6.228269936347081f, 14.17993336680509f, 56.69055260068105f};
    const float c5[3] = {4.776384997670288f, -13.74514537774601f, -65.35303263337234f};
    const float c6[3] = {-5.435455855934631f, 12.86216349749025f, 26.312873249486f};
    
    r = c0[0] + t * (c1[0] + t * (c2[0] + t * (c3[0] + t * (c4[0] + t * (c5[0] + t * c6[0])))));
    g = c0[1] + t * (c1[1] + t * (c2[1] + t * (c3[1] + t * (c4[1] + t * (c5[1] + t * c6[1])))));
    b = c0[2] + t * (c1[2] + t * (c2[2] + t * (c3[2] + t * (c4[2] + t * (c5[2] + t * c6[2])))));
    
    r = std::clamp(r, 0.0f, 1.0f);
    g = std::clamp(g, 0.0f, 1.0f);
    b = std::clamp(b, 0.0f, 1.0f);
}

void plasma(float t, float& r, float& g, float& b) {
    // Plasma colormap implementation
    const float c0[3] = {0.05873234392399702f, 0.02333670892565664f, 0.5433401826748754f};
    const float c1[3] = {2.176514634195958f, 0.2383834171260182f, 0.7539604599784036f};
    const float c2[3] = {-2.689460476458034f, -7.455851135738909f, 3.110799939717086f};
    const float c3[3] = {6.130348345893603f, 42.3461881477227f, -28.51885465332158f};
    const float c4[3] = {-11.10743619062271f, -82.66631104428044f, 60.13984767418263f};
    const float c5[3] = {10.02306557647065f, 71.41361770095349f, -54.07218655560067f};
    const float c6[3] = {-3.658713842777788f, -22.93153465461149f, 18.19190778539828f};
    
    r = c0[0] + t * (c1[0] + t * (c2[0] + t * (c3[0] + t * (c4[0] + t * (c5[0] + t * c6[0])))));
    g = c0[1] + t * (c1[1] + t * (c2[1] + t * (c3[1] + t * (c4[1] + t * (c5[1] + t * c6[1])))));
    b = c0[2] + t * (c1[2] + t * (c2[2] + t * (c3[2] + t * (c4[2] + t * (c5[2] + t * c6[2])))));
    
    r = std::clamp(r, 0.0f, 1.0f);
    g = std::clamp(g, 0.0f, 1.0f);
    b = std::clamp(b, 0.0f, 1.0f);
}

void jet(float t, float& r, float& g, float& b) {
    // Jet colormap implementation
    if (t < 0.125f) {
        r = 0.0f;
        g = 0.0f;
        b = 0.5f + 4.0f * t;
    } else if (t < 0.375f) {
        r = 0.0f;
        g = 4.0f * (t - 0.125f);
        b = 1.0f;
    } else if (t < 0.625f) {
        r = 4.0f * (t - 0.375f);
        g = 1.0f;
        b = 1.0f - 4.0f * (t - 0.375f);
    } else if (t < 0.875f) {
        r = 1.0f;
        g = 1.0f - 4.0f * (t - 0.625f);
        b = 0.0f;
    } else {
        r = 1.0f - 4.0f * (t - 0.875f);
        g = 0.0f;
        b = 0.0f;
    }
}

void grayscale(float t, float& r, float& g, float& b) {
    r = g = b = t;
}

void turbo(float t, float& r, float& g, float& b) {
    // Turbo colormap implementation
    const float c0[3] = {0.18995f, 0.07176f, 0.23217f};
    const float c1[3] = {0.30315f, 0.50427f, 0.94791f};
    const float c2[3] = {0.38343f, 0.87344f, 0.54467f};
    const float c3[3] = {0.98853f, 0.64411f, 0.09395f};
    const float c4[3] = {0.95839f, 0.82344f, 0.08425f};
    
    float t2 = t * t;
    float t3 = t2 * t;
    float t4 = t3 * t;
    
    r = c0[0] + c1[0] * t + c2[0] * t2 + c3[0] * t3 + c4[0] * t4;
    g = c0[1] + c1[1] * t + c2[1] * t2 + c3[1] * t3 + c4[1] * t4;
    b = c0[2] + c1[2] * t + c2[2] * t2 + c3[2] * t3 + c4[2] * t4;
    
    r = std::clamp(r, 0.0f, 1.0f);
    g = std::clamp(g, 0.0f, 1.0f);
    b = std::clamp(b, 0.0f, 1.0f);
}

void rdylbu(float t, float& r, float& g, float& b) {
    // Red-Yellow-Blue colormap
    if (t < 0.5f) {
        float t2 = 2.0f * t;
        r = 1.0f;
        g = t2;
        b = 0.0f;
    } else {
        float t2 = 2.0f * (t - 0.5f);
        r = 1.0f - t2;
        g = 1.0f - t2;
        b = t2;
    }
}

void apply_colormap(int colormap_id, float t, unsigned char& r, unsigned char& g, unsigned char& b) {
    float rf, gf, bf;
    
    switch (colormap_id) {
        case 0: viridis(t, rf, gf, bf); break;
        case 1: plasma(t, rf, gf, bf); break;
        case 2: jet(t, rf, gf, bf); break;
        case 3: grayscale(t, rf, gf, bf); break;
        case 4: turbo(t, rf, gf, bf); break;
        case 5: rdylbu(t, rf, gf, bf); break;
        default: viridis(t, rf, gf, bf); break;
    }
    
    r = static_cast<unsigned char>(rf * 255.0f);
    g = static_cast<unsigned char>(gf * 255.0f);
    b = static_cast<unsigned char>(bf * 255.0f);
}

void generate_colormap_texture(int colormap_id, std::vector<unsigned char>& texture_data, int width) {
    texture_data.resize(width * 3);
    for (int i = 0; i < width; ++i) {
        float t = static_cast<float>(i) / static_cast<float>(width - 1);
        apply_colormap(colormap_id, t, 
                      texture_data[i * 3], 
                      texture_data[i * 3 + 1], 
                      texture_data[i * 3 + 2]);
    }
}

} // namespace Colormaps 
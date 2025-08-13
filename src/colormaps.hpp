#pragma once

#include <vector>
#include <cstdint>

// Colormap functions for visualization
namespace Colormaps {
    
    // Convert a normalized value [0,1] to RGB color using different colormaps
    void viridis(float t, float& r, float& g, float& b);
    void plasma(float t, float& r, float& g, float& b);
    void jet(float t, float& r, float& g, float& b);
    void grayscale(float t, float& r, float& g, float& b);
    void turbo(float t, float& r, float& g, float& b);
    void rdylbu(float t, float& r, float& g, float& b);
    
    // Apply colormap to a normalized value and return RGB bytes
    void apply_colormap(int colormap_id, float t, unsigned char& r, unsigned char& g, unsigned char& b);
    
    // Generate a colormap texture for the UI
    void generate_colormap_texture(int colormap_id, std::vector<unsigned char>& texture_data, int width = 256);
}

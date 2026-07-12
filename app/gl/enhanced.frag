#version 120

varying vec2 texture_coord;
uniform sampler2D depth_tex;
uniform sampler2D shaded_tex;

uniform float zmin_local;
uniform float dz_local;
uniform float zmin_global;
uniform float dz_global;
uniform int is_2d;

uniform vec3 color;
uniform vec2 pixel_size;

// Occlusion contribution from one neighboring depth tap:
// neighbors closer to the viewer than us cast a little shadow.
float tap(vec2 d, float center)
{
    float n = texture2D(depth_tex, texture_coord + d).r;
    return clamp((n - center) * 6.0, 0.0, 1.0);
}

float occlusion(float center)
{
    vec2 r1 = pixel_size * 2.0;
    vec2 r2 = pixel_size * 5.0;
    float occ =
        tap(vec2( r1.x,  0.0), center) + tap(vec2(-r1.x,  0.0), center) +
        tap(vec2( 0.0,  r1.y), center) + tap(vec2( 0.0, -r1.y), center) +
        tap(vec2( r2.x,  r2.y), center) + tap(vec2(-r2.x,  r2.y), center) +
        tap(vec2( r2.x, -r2.y), center) + tap(vec2(-r2.x, -r2.y), center);
    return occ / 8.0;
}

vec4 shade(vec4 norm, float center_depth)
{
    vec3 light = vec3(0.99 * color.r, 0.96 * color.g, 0.89 * color.b);
    vec3 dark = vec3(0.20 * color.r, 0.25 * color.g, 0.30 * color.b);

    vec3 n = 2.0 * (norm.xyz - vec3(0.5));

    if (is_2d == 1)
    {
        float a = dot(n, vec3(0.0, 0.0, 1.0)) * 0.5 + 0.5;
        return vec4(a * light + (1.0 - a) * dark, 1.0);
    }

    // Hemispheric ambient: sky above, ground bounce below
    float hemi = n.z * 0.5 + 0.5;
    vec3 ambient = mix(dark, light, hemi) * 0.40;

    // Key light, upper-left-ish
    float key = max(dot(n, vec3(0.57, -0.57, 0.57)), 0.0);
    vec3 diffuse = key * light * 0.62;

    // Cheap screen-space AO from the depth buffer
    float ao = 1.0 - 0.5 * occlusion(center_depth);

    // Cool fresnel rim against the dark backdrop
    float rim = pow(1.0 - max(n.z, 0.0), 3.0) * 0.22;
    vec3 rim_c = rim * vec3(0.75, 0.83, 1.0);

    vec3 c = (ambient + diffuse) * ao + rim_c;

    // Gamma-correct output
    return vec4(pow(clamp(c, 0.0, 1.0), vec3(1.0 / 1.8)), 1.0);
}

void main() {
    vec4 depth = texture2D(depth_tex, texture_coord);

    if (depth.r == 0.0f)
    {
        discard;
    }

    gl_FragColor = shade(texture2D(shaded_tex, texture_coord), depth.r);

    float fd_local = depth.r * dz_local + zmin_local;
    float fd_global = (fd_local - zmin_global) / dz_global;
    gl_FragDepth = 0.9 - fd_global*0.3;
}

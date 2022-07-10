#ifndef CAMERA_H
#define CAMERA_H

#include "rtweekend.h"

class Camera
{
public:
    double fov, aspect, aperture;
    vec3 gaze, up, right;

    Camera()
    {
        auto aspect_radio = 16.0 / 9.0;
        auto viewport_height = 2.0;
        auto viewport_width = aspect_radio * viewport_height;
        auto focal_length = 1.0;

        origin = point3(0, 0, 0);
        horizontal = vec3(viewport_width, 0.0, 0.0);
        vertical = vec3(0.0, viewport_height, 0.0);
        lower_left_corner = origin - horizontal / 2 - vertical / 2 - vec3(0, 0, focal_length);
    }

    Camera(double f, double a): fov(f), aspect(a) 
    {
        double theta = degree_to_radius(fov);
        double focal_length = 1.0;
        double viewport_height = 2.0 * tan(theta / 2.0) * focal_length;
        double viewport_width = aspect * viewport_height;

        origin = point3(0, 0, 0);
        horizontal = vec3(viewport_width, 0.0, 0.0);
        vertical = vec3(0.0, viewport_height, 0.0);
        lower_left_corner = origin - horizontal / 2 - vertical / 2 - vec3(0, 0, focal_length);
    }

    Camera(double f, double a, point3 camera_pos, vec3 up): fov(f), aspect(a), camera_position(camera_pos)
    {
        double theta = degree_to_radius(fov);
        double camera_origin_dist = (origin - camera_pos).length();
        double viewport_height = 2.0 * tan(theta / 2.0) * camera_origin_dist;
        double viewport_width = aspect * viewport_height;

        origin = point3(0, 0, 0);
        g_ = unit_vector(origin - camera_pos);
        v_ = unit_vector(cross(g_, up));
        u_ = unit_vector(cross(v_, g_));

        horizontal = viewport_width * v_;
        vertical = viewport_height * u_;
        lower_left_corner = origin - horizontal / 2 - vertical / 2;
    }

    Camera(double f, double asp, double ape, point3 camera_pos, vec3 up): fov(f), aspect(asp), aperture(ape), camera_position(camera_pos)
    {
        double theta = degree_to_radius(fov);
        double camera_origin_dist = (origin - camera_pos).length();
        double viewport_height = 2.0 * tan(theta / 2.0) * camera_origin_dist;
        double viewport_width = aspect * viewport_height;
        len_radius = aperture / 2.0;

        origin = point3(0, 0, 0);
        g_ = unit_vector(origin - camera_pos);
        v_ = unit_vector(cross(g_, up));
        u_ = unit_vector(cross(v_, g_));

        horizontal = viewport_width * v_;
        vertical = viewport_height * u_;
        lower_left_corner = origin - horizontal / 2 - vertical / 2;
    }

    ray get_ray(double u, double v) const
    {
        if(camera_position == origin)
        {
            return ray(origin, lower_left_corner + u * horizontal + v * vertical - origin);
        }
        else
        {
            vec3 rd = len_radius * random_in_unit_disk();
            vec3 offset = v_ * rd.x() + u_ * rd.y();
            return ray(camera_position + offset, lower_left_corner + u * horizontal + v * vertical - camera_position - offset);
        }
    }

private:
    point3 origin;
    point3 lower_left_corner;
    point3 camera_position;
    vec3 horizontal;
    vec3 vertical;
    vec3 u_, v_, g_;
    double len_radius = 0.0;
};


#endif
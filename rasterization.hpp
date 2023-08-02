#ifndef RASTERIZATION_HPP
#define RASTERIZATION_HPP

#include <cmath>
#include <vector>
#include <array>
#include <functional>
#include "thread_pool_manager.hpp"
#include "bitmap.h"

using namespace bitmap;

namespace pipeline3D
{

    template <class VertexData>
    struct Vertex
    {
        double x;
        double y;
        double z;
        VertexData vertexData;

        Vertex() = default;
        Vertex(const double &_x, const double &_y, const double &_z, const VertexData &_data) : x(_x), y(_y), z(_z), vertexData(_data) {}
    };

    /**
     * @brief Triangle 3D object, which will wrap its 3 vertices.
     *
     */
    template <class VertexData>
    class Triangle3D
    {
    public:
        /**
         * @brief Triangle's vertices.
         *
         */
        Vertex<VertexData> v0, v1, v2;

        /**
         * @brief Default constructor
         *
         */
        Triangle3D() = default;

        /**
         * @brief Construct a new Triangle 3D object, whose vertices will be those passed as parameters.
         *
         * @param v0
         * @param v1
         * @param v2
         */
        Triangle3D(Vertex<VertexData> _v0, Vertex<VertexData> _v1, Vertex<VertexData> _v2) : v0(_v0), v1(_v1), v2(_v2) {}

        Triangle3D(Vertex<VertexData> &&_v0, Vertex<VertexData> &&_v1, Vertex<VertexData> &&_v2) : v0(_v0), v1(_v1), v2(_v2)
        {
            _v0 = 0;
            _v1 = 0;
            _v2 = 0;
        }

        /**
         * @brief Default destructor.
         *
         */
        ~Triangle3D() = default;
    };

    template <class Target_t, class VertexData, class Derived>
    struct VertexShader
    {
        Target_t shade(Vertex<VertexData> v)
        {
            static_assert(std::is_base_of<VertexShader, Derived>::value,
                          "Implemented shader must inherit from VertexShader");
            return static_cast<Derived *>(this)->shade(v);
        }

        Vertex<VertexData> interpolate(const Vertex<VertexData> &v1,
                                       const Vertex<VertexData> &v2,
                                       double w)
        {
            static_assert(std::is_base_of<VertexShader, Derived>::value,
                          "Implemented shader must inherit from VertexShader");
            return static_cast<Derived *>(this)->interpolate(v1, v2, w);
        }

        void perspectiveCorrect(Vertex<VertexData> &v)
        {
            static_assert(std::is_base_of<VertexShader, Derived>::value,
                          "Implemented shader must inherit from VertexShader");
            return static_cast<Derived *>(this)->perspectiveCorrect(v);
        }

        ~VertexShader() = default;
    };

    template <class Target_t, class ZBuffer_t>
    class Rasterizer
    {
    public:
        void set_target(int w, int h, Target_t *t, ZBuffer_t *z)
        {
            width = w;
            height = h;
            target = t;
            z_buffer = z;
        }

        std::vector<Target_t> get_z_buffer() { return std::move(z_buffer); }

        void set_perspective_projection(double left, double right, double top, double bottom, double near, double far)
        {
            const double w = right - left;
            const double h = bottom - top;
            const double d = far - near;

            // row-major
            projection_matrix[0] = 2.0f * near / w;
            projection_matrix[1] = 0;
            projection_matrix[2] = -(right + left) / w;
            projection_matrix[3] = 0;
            projection_matrix[4 * 1 + 0] = 0;
            projection_matrix[4 * 1 + 1] = 2.0f * near / h;
            projection_matrix[4 * 1 + 2] = -(bottom + top) / h;
            projection_matrix[4 * 1 + 3] = 0;
            projection_matrix[4 * 2 + 0] = 0;
            projection_matrix[4 * 2 + 1] = 0;
            projection_matrix[4 * 2 + 2] = (far + near) / d;
            projection_matrix[4 * 2 + 3] = -2.0f * far * near / d;
            projection_matrix[4 * 3 + 0] = 0;
            projection_matrix[4 * 3 + 1] = 0;
            projection_matrix[4 * 3 + 2] = 1;
            projection_matrix[4 * 3 + 3] = 0;
        }

        void set_orthographic_projection(double left, double right, double top, double bottom, double near, double far)
        {
            const double w = right - left;
            const double h = bottom - top;
            const double d = far - near;

            // row-major
            projection_matrix[0] = 2.0f / w;
            projection_matrix[1] = 0;
            projection_matrix[2] = 0;
            projection_matrix[3] = -(right + left) / w;
            projection_matrix[4 * 1 + 0] = 0;
            projection_matrix[4 * 1 + 1] = 2.0f / h;
            projection_matrix[4 * 1 + 2] = 0;
            projection_matrix[4 * 1 + 3] = -(bottom + top) / h;
            projection_matrix[4 * 2 + 0] = 0;
            projection_matrix[4 * 2 + 1] = 0;
            projection_matrix[4 * 2 + 2] = 2.0f / d;
            projection_matrix[4 * 2 + 3] = -(far + near) / d;
            projection_matrix[4 * 3 + 0] = 0;
            projection_matrix[4 * 3 + 1] = 0;
            projection_matrix[4 * 3 + 2] = 0;
            projection_matrix[4 * 3 + 3] = 1;
        }

        template <class Shader, class VertexData>
        void render_triangle(Vertex<VertexData> v1,
                             Vertex<VertexData> v2,
                             Vertex<VertexData> v3,
                             Shader &shader)
        {

            // project view coordinates to ndc;
            std::array<double, 3> ndc_v1;
            std::array<double, 3> ndc_v2;
            std::array<double, 3> ndc_v3;
            project(v1, ndc_v1);
            project(v2, ndc_v2);
            project(v3, ndc_v3);
            // at first sort the three vertices by y-coordinate ascending so v1 is the topmost vertice
            if (ndc_v1[1] > ndc_v2[1])
            {
                std::swap(v1, v2);
                std::swap(ndc_v1, ndc_v2);
            }
            if (ndc_v1[1] > ndc_v3[1])
            {
                std::swap(v1, v3);
                std::swap(ndc_v1, ndc_v3);
            }
            if (ndc_v2[1] > ndc_v3[1])
            {
                std::swap(v2, v3);
                std::swap(ndc_v2, ndc_v3);
            }
            shader.perspectiveCorrect(v1);
            shader.perspectiveCorrect(v2);
            shader.perspectiveCorrect(v3);
            // convert normalized device coordinates into Pixel coordinates
            const double x1f = ndc2idxf(ndc_v1[0], width);
            const int x1 = static_cast<int>(x1f + 0.5f);
            const double y1f = ndc2idxf(ndc_v1[1], height);
            const int y1 = static_cast<int>(y1f);
            const double x2f = ndc2idxf(ndc_v2[0], width);
            const int x2 = static_cast<int>(x2f);
            const double y2f = ndc2idxf(ndc_v2[1], height);
            const int y2 = static_cast<int>(y2f + 0.5f);
            const double x3f = ndc2idxf(ndc_v3[0], width);
            const int x3 = static_cast<int>(x3f + 0.5f);
            const double y3f = ndc2idxf(ndc_v3[1], height);
            const int y3 = static_cast<int>(y3f);

            const double idy12 = 1.0f / (y2f - y1f);
            const double idy13 = 1.0f / (y3f - y1f);
            const double idy23 = 1.0f / (y3f - y2f);

            const double m12 = (x2f - x1f) * idy12;
            const double m13 = (x3f - x1f) * idy13;
            const double m23 = (x3f - x2f) * idy23;

            const double q12 = (x1f * y2f - x2f * y1f) * idy12;
            const double q13 = (x1f * y3f - x3f * y1f) * idy13;
            const double q23 = (x2f * y3f - x3f * y2f) * idy23;

            bool horizontal12 = std::abs(m12) > 1.0f;
            bool horizontal13 = std::abs(m13) > 1.0f;
            bool horizontal23 = std::abs(m23) > 1.0f;

            if (y1 >= height || y3 < 0)
                return;

            if (m13 > m12)
            { // v2 is on the left of the line v1-v3
                int y = y1;
                double w1f = (y2f - y) * idy12;
                double w1l = (y3f - y) * idy13;
                double xf = m12 * y + q12;
                double xl = m13 * y + q13;
                while (y != y2)
                {
                    const int first = horizontal12 ? static_cast<int>(xf + 0.5 * m12) : static_cast<int>(xf);
                    const int last = horizontal13 ? static_cast<int>(xl + 0.5 * m13) + 1 : static_cast<int>(xl) + 1;

                    const double step = 1.0f / (xl - xf);
                    const double w0 = 1.0f + (xf - first) * step;
                    Vertex<VertexData> vl = shader.interpolate(v1, v2, w1f);
                    Vertex<VertexData> vr = shader.interpolate(v1, v3, w1l);

                    render_scanline<Shader, VertexData>(y,
                                                        first,
                                                        last,
                                                        vl,
                                                        vr,
                                                        interpolateF(ndc_v1[2], ndc_v2[2], w1f),
                                                        interpolateF(ndc_v1[2], ndc_v3[2], w1l),
                                                        w0,
                                                        step,
                                                        shader);
                    ++y;
                    if (y == height)
                    {
                        return;
                    }
                    w1f -= idy12;
                    w1l -= idy13;
                    xf += m12;
                    xl += m13;
                }
                if (y2 == y3)
                {
                    const double step = 1.0f / (xl - xf);
                    const double w0 = 1.0f + (xf - x2) * step;
                    render_scanline<Shader, VertexData>(y,
                                                        x2,
                                                        x3,
                                                        v2,
                                                        v3,
                                                        ndc_v2[2],
                                                        ndc_v3[2],
                                                        w0,
                                                        step,
                                                        shader);

                    return;
                }
                else
                {
                    if (std::abs(m12) > std::abs(m23))
                    {
                        xf = m23 * y + q23;
                    }
                    const double step = 1.0f / (xl - xf);
                    const double w0 = 1.0f + (xf - x2) * step;
                    const int last = horizontal13 ? static_cast<int>(xl + 0.5 * m13) + 1 : static_cast<int>(xl) + 1;

                    render_scanline<Shader, VertexData>(y,
                                                        x2,
                                                        last,
                                                        v2,
                                                        shader.interpolate(v1, v3, w1l),
                                                        ndc_v2[2],
                                                        interpolateF(ndc_v1[2], ndc_v3[2], w1l),
                                                        w0,
                                                        step,
                                                        shader);
                }
                ++y;
                if (y == height)
                {
                    return;
                }
                w1l -= idy13;
                double w2f = (y3f - y) * idy23;
                xf = m23 * y + q23;
                xl += m13;

                while (y != y3)
                {
                    const int first = horizontal23 ? static_cast<int>(xf + 0.5 * m23) : static_cast<int>(xf);
                    const int last = horizontal13 ? static_cast<int>(xl + 0.5 * m13) + 1 : static_cast<int>(xl) + 1;

                    const double step = 1.0f / (xl - xf);
                    const double w0 = 1.0f + (xf - first) * step;

                    render_scanline<Shader, VertexData>(y,
                                                        first,
                                                        last,
                                                        shader.interpolate(v2, v3, w2f),
                                                        shader.interpolate(v1, v3, w1l),
                                                        interpolateF(ndc_v2[2], ndc_v3[2], w2f),
                                                        interpolateF(ndc_v1[2], ndc_v3[2], w1l),
                                                        w0,
                                                        step,
                                                        shader);
                    ++y;
                    if (y == height)
                    {
                        return;
                    }
                    w1l -= idy13;
                    w2f -= idy23;
                    xf += m23;
                    xl += m13;
                }

                const int first = horizontal23 ? static_cast<int>(xf + 0.5 * m23) : static_cast<int>(xf);
                const double step = 1.0f / (xl - xf);
                const double w0 = 1.0f + (xf - first) * step;
                render_scanline<Shader, VertexData>(y,
                                                    first,
                                                    x3,
                                                    shader.interpolate(v2, v3, w2f), v3,
                                                    interpolateF(ndc_v2[2], ndc_v3[2], w2f),
                                                    ndc_v3[2],
                                                    w0,
                                                    step,
                                                    shader);
            }
            else
            { // v2 is on the right of the line v1-v3
                int y = y1;
                double w1f = (y3f - y) * idy13;
                double w1l = (y2f - y) * idy12;
                double xf = m13 * y + q13;
                double xl = m12 * y + q12;
                while (y != y2)
                {
                    const int first = horizontal13 ? static_cast<int>(xf + 0.5 * m13) : static_cast<int>(xf);
                    const int last = horizontal12 ? static_cast<int>(xl + 0.5 * m12) + 1 : static_cast<int>(xl) + 1;

                    const double step = 1.0f / (xl - xf);
                    const double w0 = 1.0f + (xf - first) * step;

                    render_scanline<Shader, VertexData>(y,
                                                        first,
                                                        last,
                                                        shader.interpolate(v1, v3, w1f),
                                                        shader.interpolate(v1, v2, w1l),
                                                        interpolateF(ndc_v1[2], ndc_v3[2], w1f),
                                                        interpolateF(ndc_v1[2], ndc_v2[2], w1l),
                                                        w0,
                                                        step,
                                                        shader);
                    ++y;
                    if (y == height)
                    {
                        return;
                    }
                    w1f -= idy13;
                    w1l -= idy12;
                    xf += m13;
                    xl += m12;
                }
                if (y2 == y3)
                {
                    const double step = 1.0f / (xl - xf);
                    const double w0 = 1.0f + (xf - x3) * step;
                    render_scanline<Shader, VertexData>(y,
                                                        x3,
                                                        x2 + 1,
                                                        v3,
                                                        v2,
                                                        ndc_v3[2],
                                                        ndc_v2[2],
                                                        w0,
                                                        step,
                                                        shader);

                    return;
                }
                else
                {
                    const int first = horizontal13 ? static_cast<int>(xf + 0.5 * m13) : static_cast<int>(xf);
                    const double step = 1.0f / (x2 + 1 - xf);
                    const double w0 = 1.0f + (xf - first) * step;

                    render_scanline<Shader, VertexData>(y,
                                                        first,
                                                        x2 + 1,
                                                        shader.interpolate(v1, v3, w1f),
                                                        v2,
                                                        interpolateF(ndc_v1[2], ndc_v3[2], w1f),
                                                        ndc_v2[2],
                                                        w0,
                                                        step,
                                                        shader);
                }
                ++y;
                if (y == height)
                {
                    return;
                }
                w1f -= idy13;
                double w2l = (y3f - y) * idy23;
                xl = m23 * y + q23;
                xf += m13;

                while (y != y3)
                {
                    const int first = horizontal13 ? static_cast<int>(xf + 0.5 * m13) : static_cast<int>(xf);
                    const int last = horizontal23 ? static_cast<int>(xl + 0.5 * m23) + 1 : static_cast<int>(xl) + 1;

                    const double step = 1.0f / (xl - xf);
                    const double w0 = 1.0f + (xf - first) * step;

                    render_scanline<Shader, VertexData>(y, //-> 10
                                                        first,
                                                        last,
                                                        shader.interpolate(v1, v3, w1f),
                                                        shader.interpolate(v2, v3, w2l),
                                                        interpolateF(ndc_v1[2], ndc_v3[2], w1f),
                                                        interpolateF(ndc_v2[2], ndc_v3[2], w2l),
                                                        w0,
                                                        step,
                                                        shader);
                    ++y;
                    if (y == height)
                    {
                        return;
                    }
                    w1f -= idy13;
                    w2l -= idy23;
                    xl += m23;
                    xf += m13;
                }
                const double step = 1.0f / (xl - xf);
                const double w0 = 1.0f + (xf - x3) * step;
                const int last = horizontal23 ? static_cast<int>(xl + 0.5 * m23) + 1 : static_cast<int>(xl) + 1;

                render_scanline<Shader, VertexData>(y,
                                                    x3,
                                                    last,
                                                    v3,
                                                    shader.interpolate(v2, v3, w2l),
                                                    ndc_v3[2],
                                                    interpolateF(ndc_v2[2], ndc_v3[2], w2l),
                                                    w0,
                                                    step,
                                                    shader);
            }
        }

        template <class Shader, class VertexData>
        void render_triangle(Vertex<VertexData> v1,
                             Vertex<VertexData> v2,
                             Vertex<VertexData> v3,
                             Shader &shader,
                             thread_pool_manager<std::mutex *> &threadPool,
                             const bool &isScanlineMultithreaded)
        {
            // project view coordinates to ndc;
            std::array<double, 3> ndc_v1;
            std::array<double, 3> ndc_v2;
            std::array<double, 3> ndc_v3;
            project(v1, ndc_v1);
            project(v2, ndc_v2);
            project(v3, ndc_v3);
            // at first sort the three vertices by y-coordinate ascending so v1 is the topmost vertice
            if (ndc_v1[1] > ndc_v2[1])
            {
                std::swap(v1, v2);
                std::swap(ndc_v1, ndc_v2);
            }
            if (ndc_v1[1] > ndc_v3[1])
            {
                std::swap(v1, v3);
                std::swap(ndc_v1, ndc_v3);
            }
            if (ndc_v2[1] > ndc_v3[1])
            {
                std::swap(v2, v3);
                std::swap(ndc_v2, ndc_v3);
            }
            shader.perspectiveCorrect(v1);
            shader.perspectiveCorrect(v2);
            shader.perspectiveCorrect(v3);

            // convert normalized device coordinates into Pixel coordinates
            const double x1f = ndc2idxf(ndc_v1[0], width);
            const int x1 = static_cast<int>(x1f + 0.5f);
            const double y1f = ndc2idxf(ndc_v1[1], height);
            const int y1 = static_cast<int>(y1f);
            const double x2f = ndc2idxf(ndc_v2[0], width);
            const int x2 = static_cast<int>(x2f);
            const double y2f = ndc2idxf(ndc_v2[1], height);
            const int y2 = static_cast<int>(y2f + 0.5f);
            const double x3f = ndc2idxf(ndc_v3[0], width);
            const int x3 = static_cast<int>(x3f + 0.5f);
            const double y3f = ndc2idxf(ndc_v3[1], height);
            const int y3 = static_cast<int>(y3f);

            const double idy12 = 1.0f / (y2f - y1f);
            const double idy13 = 1.0f / (y3f - y1f);
            const double idy23 = 1.0f / (y3f - y2f);

            const double m12 = (x2f - x1f) * idy12;
            const double m13 = (x3f - x1f) * idy13;
            const double m23 = (x3f - x2f) * idy23;

            const double q12 = (x1f * y2f - x2f * y1f) * idy12;
            const double q13 = (x1f * y3f - x3f * y1f) * idy13;
            const double q23 = (x2f * y3f - x3f * y2f) * idy23;

            bool horizontal12 = std::abs(m12) > 1.0f;
            bool horizontal13 = std::abs(m13) > 1.0f;
            bool horizontal23 = std::abs(m23) > 1.0f;

            int first, last;
            double step, w0;

            if (isScanlineMultithreaded)
            {
                if (m13 > m12)
                { // v2 is on the left of the line v1-v3
                    int y = y1;
                    double w1f = (y2f - y) * idy12;
                    double w1l = (y3f - y) * idy13;
                    double xf = m12 * y + q12;
                    double xl = m13 * y + q13;

                    while (y != y2)
                    {
                        first = horizontal12 ? static_cast<int>(xf + 0.5 * m12) : static_cast<int>(xf);
                        last = horizontal13 ? static_cast<int>(xl + 0.5 * m13) + 1 : static_cast<int>(xl) + 1;

                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - first) * step;
                        Vertex<VertexData> vl = shader.interpolate(v1, v2, w1f);
                        Vertex<VertexData> vr = shader.interpolate(v1, v3, w1l);

                        double interp1, interp2;
                        interp1 = interpolateF(ndc_v1[2], ndc_v2[2], w1f);
                        interp2 = interpolateF(ndc_v1[2], ndc_v3[2], w1l);

                        threadPool.executeTask(std::forward<std::function<void()>>([&,
                                                                                    y,
                                                                                    first,
                                                                                    last,
                                                                                    vl,
                                                                                    vr,
                                                                                    interp1,
                                                                                    interp2,
                                                                                    w0,
                                                                                    step]
                                                                                   { render_scanline<Shader, VertexData>(y,
                                                                                                                         first,
                                                                                                                         last,
                                                                                                                         vl,
                                                                                                                         vr,
                                                                                                                         interp1,
                                                                                                                         interp2,
                                                                                                                         w0,
                                                                                                                         step,
                                                                                                                         shader,
                                                                                                                         threadPool); }));
                        ++y;
                        if (y == height)
                        {
                            return;
                        }
                        w1f -= idy12;
                        w1l -= idy13;
                        xf += m12;
                        xl += m13;
                    }
                    if (y2 == y3)
                    {
                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - x2) * step;

                        threadPool.executeTask(std::forward<std::function<void()>>([&,
                                                                                    y,
                                                                                    x2,
                                                                                    x3,
                                                                                    v2,
                                                                                    v3,
                                                                                    ndc_v2,
                                                                                    ndc_v3,
                                                                                    w0,
                                                                                    step]
                                                                                   { render_scanline<Shader, VertexData>(y,
                                                                                                                         x2,
                                                                                                                         x3,
                                                                                                                         v2,
                                                                                                                         v3,
                                                                                                                         ndc_v2[2],
                                                                                                                         ndc_v3[2],
                                                                                                                         w0,
                                                                                                                         step,
                                                                                                                         shader,
                                                                                                                         threadPool); }));
                        return;
                    }
                    else
                    {
                        if (std::abs(m12) > std::abs(m23))
                        {
                            xf = m23 * y + q23;
                        }

                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - x2) * step;
                        last = horizontal13 ? static_cast<int>(xl + 0.5 * m13) + 1 : static_cast<int>(xl) + 1;

                        threadPool.executeTask(std::forward<std::function<void()>>([&,
                                                                                    y,
                                                                                    last,
                                                                                    x2,
                                                                                    v1,
                                                                                    v2,
                                                                                    v3,
                                                                                    w1l,
                                                                                    ndc_v1,
                                                                                    ndc_v2,
                                                                                    ndc_v3,
                                                                                    w0,
                                                                                    step]
                                                                                   { render_scanline<Shader, VertexData>(y,
                                                                                                                         x2,
                                                                                                                         last,
                                                                                                                         v2,
                                                                                                                         shader.interpolate(v1, v3, w1l),
                                                                                                                         ndc_v2[2],
                                                                                                                         interpolateF(ndc_v1[2], ndc_v3[2], w1l),
                                                                                                                         w0,
                                                                                                                         step,
                                                                                                                         shader,
                                                                                                                         threadPool); }));
                    }
                    ++y;
                    if (y == height)
                    {
                        return;
                    }
                    w1l -= idy13;
                    double w2f = (y3f - y) * idy23;
                    xf = m23 * y + q23;
                    xl += m13;

                    while (y != y3)
                    {
                        first = horizontal23 ? static_cast<int>(xf + 0.5 * m23) : static_cast<int>(xf);
                        last = horizontal13 ? static_cast<int>(xl + 0.5 * m13) + 1 : static_cast<int>(xl) + 1;

                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - first) * step;

                        threadPool.executeTask(std::forward<std::function<void()>>([&,
                                                                                    y,
                                                                                    first,
                                                                                    last,
                                                                                    v1,
                                                                                    v2,
                                                                                    v3,
                                                                                    w1l,
                                                                                    w2f,
                                                                                    ndc_v1,
                                                                                    ndc_v2,
                                                                                    ndc_v3,
                                                                                    w0,
                                                                                    step]
                                                                                   { render_scanline<Shader, VertexData>(y,
                                                                                                                         first,
                                                                                                                         last,
                                                                                                                         shader.interpolate(v2, v3, w2f),
                                                                                                                         shader.interpolate(v1, v3, w1l),
                                                                                                                         interpolateF(ndc_v2[2], ndc_v3[2], w2f),
                                                                                                                         interpolateF(ndc_v1[2], ndc_v3[2], w1l),
                                                                                                                         w0,
                                                                                                                         step,
                                                                                                                         shader,
                                                                                                                         threadPool); }));
                        ++y;
                        if (y == height)
                        {
                            return;
                        }
                        w1l -= idy13;
                        w2f -= idy23;
                        xf += m23;
                        xl += m13;
                    }

                    first = horizontal23 ? static_cast<int>(xf + 0.5 * m23) : static_cast<int>(xf);
                    step = 1.0f / (xl - xf);
                    w0 = 1.0f + (xf - first) * step;
                    threadPool.executeTask(std::forward<std::function<void()>>([&,
                                                                                y,
                                                                                first,
                                                                                x3,
                                                                                v2,
                                                                                v3,
                                                                                w2f,
                                                                                ndc_v2,
                                                                                ndc_v3,
                                                                                w0,
                                                                                step]
                                                                               { render_scanline<Shader, VertexData>(y,
                                                                                                                     first,
                                                                                                                     x3,
                                                                                                                     shader.interpolate(v2, v3, w2f), v3,
                                                                                                                     interpolateF(ndc_v2[2], ndc_v3[2], w2f),
                                                                                                                     ndc_v3[2],
                                                                                                                     w0,
                                                                                                                     step,
                                                                                                                     shader,
                                                                                                                     threadPool); }));
                }
                else
                { // v2 is on the right of the line v1-v3
                    int y = y1;
                    double w1f = (y3f - y) * idy13;
                    double w1l = (y2f - y) * idy12;
                    double xf = m13 * y + q13;
                    double xl = m12 * y + q12;
                    while (y != y2)
                    {
                        first = horizontal13 ? static_cast<int>(xf + 0.5 * m13) : static_cast<int>(xf);
                        last = horizontal12 ? static_cast<int>(xl + 0.5 * m12) + 1 : static_cast<int>(xl) + 1;

                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - first) * step;

                        threadPool.executeTask(std::forward<std::function<void()>>([&,
                                                                                    y,
                                                                                    first,
                                                                                    last,
                                                                                    v1,
                                                                                    v2,
                                                                                    v3,
                                                                                    w1l,
                                                                                    w1f,
                                                                                    ndc_v1,
                                                                                    ndc_v2,
                                                                                    ndc_v3,
                                                                                    w0,
                                                                                    step]
                                                                                   { render_scanline<Shader, VertexData>(y,
                                                                                                                         first,
                                                                                                                         last,
                                                                                                                         shader.interpolate(v1, v3, w1f),
                                                                                                                         shader.interpolate(v1, v2, w1l),
                                                                                                                         interpolateF(ndc_v1[2], ndc_v3[2], w1f),
                                                                                                                         interpolateF(ndc_v1[2], ndc_v2[2], w1l),
                                                                                                                         w0,
                                                                                                                         step,
                                                                                                                         shader,
                                                                                                                         threadPool); }));
                        ++y;
                        if (y == height)
                        {
                            return;
                        }
                        w1f -= idy13;
                        w1l -= idy12;
                        xf += m13;
                        xl += m12;
                    }
                    if (y2 == y3)
                    {
                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - x3) * step;
                        threadPool.executeTask(std::forward<std::function<void()>>([&,
                                                                                    y,
                                                                                    x2,
                                                                                    x3,
                                                                                    v2,
                                                                                    v3,
                                                                                    ndc_v2,
                                                                                    ndc_v3,
                                                                                    w0,
                                                                                    step]
                                                                                   { render_scanline<Shader, VertexData>(y,
                                                                                                                         x3,
                                                                                                                         x2 + 1,
                                                                                                                         v3,
                                                                                                                         v2,
                                                                                                                         ndc_v3[2],
                                                                                                                         ndc_v2[2],
                                                                                                                         w0,
                                                                                                                         step,
                                                                                                                         shader,
                                                                                                                         threadPool); }));
                        return;
                    }
                    else
                    {
                        first = horizontal13 ? static_cast<int>(xf + 0.5 * m13) : static_cast<int>(xf);
                        step = 1.0f / (x2 + 1 - xf);
                        w0 = 1.0f + (xf - first) * step;

                        threadPool.executeTask(std::forward<std::function<void()>>([&,
                                                                                    y,
                                                                                    first,
                                                                                    x2,
                                                                                    v1,
                                                                                    v2,
                                                                                    v3,
                                                                                    w1f,
                                                                                    ndc_v1,
                                                                                    ndc_v2,
                                                                                    ndc_v3,
                                                                                    w0,
                                                                                    step]
                                                                                   { render_scanline<Shader, VertexData>(y,
                                                                                                                         first,
                                                                                                                         x2 + 1,
                                                                                                                         shader.interpolate(v1, v3, w1f),
                                                                                                                         v2,
                                                                                                                         interpolateF(ndc_v1[2], ndc_v3[2], w1f),
                                                                                                                         ndc_v2[2],
                                                                                                                         w0,
                                                                                                                         step,
                                                                                                                         shader,
                                                                                                                         threadPool); }));
                    }
                    ++y;
                    if (y == height)
                    {
                        return;
                    }
                    w1f -= idy13;
                    double w2l = (y3f - y) * idy23;
                    xl = m23 * y + q23;
                    xf += m13;

                    while (y != y3)
                    {
                        first = horizontal13 ? static_cast<int>(xf + 0.5 * m13) : static_cast<int>(xf);
                        last = horizontal23 ? static_cast<int>(xl + 0.5 * m23) + 1 : static_cast<int>(xl) + 1;

                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - first) * step;

                        threadPool.executeTask(std::forward<std::function<void()>>([&,
                                                                                    y,
                                                                                    first,
                                                                                    last,
                                                                                    v1,
                                                                                    v2,
                                                                                    v3,
                                                                                    w1f,
                                                                                    w2l,
                                                                                    ndc_v1,
                                                                                    ndc_v2,
                                                                                    ndc_v3,
                                                                                    w0,
                                                                                    step]
                                                                                   { render_scanline<Shader, VertexData>(y, //-> 10
                                                                                                                         first,
                                                                                                                         last,
                                                                                                                         shader.interpolate(v1, v3, w1f),
                                                                                                                         shader.interpolate(v2, v3, w2l),
                                                                                                                         interpolateF(ndc_v1[2], ndc_v3[2], w1f),
                                                                                                                         interpolateF(ndc_v2[2], ndc_v3[2], w2l),
                                                                                                                         w0,
                                                                                                                         step,
                                                                                                                         shader,
                                                                                                                         threadPool); }));
                        ++y; // -> 11
                        if (y == height)
                        {
                            return;
                        }
                        w1f -= idy13;
                        w2l -= idy23;
                        xl += m23;
                        xf += m13;
                    }
                    step = 1.0f / (xl - xf);
                    w0 = 1.0f + (xf - x3) * step;
                    last = horizontal23 ? static_cast<int>(xl + 0.5 * m23) + 1 : static_cast<int>(xl) + 1;

                    threadPool.executeTask(std::forward<std::function<void()>>([&,
                                                                                y,
                                                                                last,
                                                                                x3,
                                                                                v2,
                                                                                v3,
                                                                                w2l,
                                                                                ndc_v2,
                                                                                ndc_v3,
                                                                                w0,
                                                                                step]
                                                                               { render_scanline<Shader, VertexData>(y,
                                                                                                                     x3,
                                                                                                                     last,
                                                                                                                     v3,
                                                                                                                     shader.interpolate(v2, v3, w2l),
                                                                                                                     ndc_v3[2],
                                                                                                                     interpolateF(ndc_v2[2], ndc_v3[2], w2l),
                                                                                                                     w0,
                                                                                                                     step,
                                                                                                                     shader,
                                                                                                                     threadPool); }));
                }
            }
            else
            {
                if (m13 > m12)
                { // v2 is on the left of the line v1-v3
                    int y = y1;
                    double w1f = (y2f - y) * idy12;
                    double w1l = (y3f - y) * idy13;
                    double xf = m12 * y + q12;
                    double xl = m13 * y + q13;

                    while (y != y2)
                    {
                        first = horizontal12 ? static_cast<int>(xf + 0.5 * m12) : static_cast<int>(xf);
                        last = horizontal13 ? static_cast<int>(xl + 0.5 * m13) + 1 : static_cast<int>(xl) + 1;

                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - first) * step;
                        Vertex<VertexData> vl = shader.interpolate(v1, v2, w1f);
                        Vertex<VertexData> vr = shader.interpolate(v1, v3, w1l);

                        double interp1, interp2;
                        interp1 = interpolateF(ndc_v1[2], ndc_v2[2], w1f);
                        interp2 = interpolateF(ndc_v1[2], ndc_v3[2], w1l);

                        render_scanline<Shader, VertexData>(y,
                                                            first,
                                                            last,
                                                            vl,
                                                            vr,
                                                            interp1,
                                                            interp2,
                                                            w0,
                                                            step,
                                                            shader,
                                                            threadPool);
                        ++y;
                        if (y == height)
                        {
                            return;
                        }
                        w1f -= idy12;
                        w1l -= idy13;
                        xf += m12;
                        xl += m13;
                    }
                    if (y2 == y3)
                    {
                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - x2) * step;

                        render_scanline<Shader, VertexData>(y,
                                                            x2,
                                                            x3,
                                                            v2,
                                                            v3,
                                                            ndc_v2[2],
                                                            ndc_v3[2],
                                                            w0,
                                                            step,
                                                            shader,
                                                            threadPool);
                        return;
                    }
                    else
                    {
                        if (std::abs(m12) > std::abs(m23))
                        {
                            xf = m23 * y + q23;
                        }

                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - x2) * step;
                        last = horizontal13 ? static_cast<int>(xl + 0.5 * m13) + 1 : static_cast<int>(xl) + 1;

                        render_scanline<Shader, VertexData>(y,
                                                            x2,
                                                            last,
                                                            v2,
                                                            shader.interpolate(v1, v3, w1l),
                                                            ndc_v2[2],
                                                            interpolateF(ndc_v1[2], ndc_v3[2], w1l),
                                                            w0,
                                                            step,
                                                            shader,
                                                            threadPool);
                    }
                    ++y;
                    if (y == height)
                    {
                        return;
                    }
                    w1l -= idy13;
                    double w2f = (y3f - y) * idy23;
                    xf = m23 * y + q23;
                    xl += m13;

                    while (y != y3)
                    {
                        first = horizontal23 ? static_cast<int>(xf + 0.5 * m23) : static_cast<int>(xf);
                        last = horizontal13 ? static_cast<int>(xl + 0.5 * m13) + 1 : static_cast<int>(xl) + 1;

                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - first) * step;

                        render_scanline<Shader, VertexData>(y,
                                                            first,
                                                            last,
                                                            shader.interpolate(v2, v3, w2f),
                                                            shader.interpolate(v1, v3, w1l),
                                                            interpolateF(ndc_v2[2], ndc_v3[2], w2f),
                                                            interpolateF(ndc_v1[2], ndc_v3[2], w1l),
                                                            w0,
                                                            step,
                                                            shader,
                                                            threadPool);
                        ++y;
                        if (y == height)
                        {
                            return;
                        }
                        w1l -= idy13;
                        w2f -= idy23;
                        xf += m23;
                        xl += m13;
                    }

                    first = horizontal23 ? static_cast<int>(xf + 0.5 * m23) : static_cast<int>(xf);
                    step = 1.0f / (xl - xf);
                    w0 = 1.0f + (xf - first) * step;
                    render_scanline<Shader, VertexData>(y,
                                                        first,
                                                        x3,
                                                        shader.interpolate(v2, v3, w2f), v3,
                                                        interpolateF(ndc_v2[2], ndc_v3[2], w2f),
                                                        ndc_v3[2],
                                                        w0,
                                                        step,
                                                        shader,
                                                        threadPool);
                }
                else
                { // v2 is on the right of the line v1-v3
                    int y = y1;
                    double w1f = (y3f - y) * idy13;
                    double w1l = (y2f - y) * idy12;
                    double xf = m13 * y + q13;
                    double xl = m12 * y + q12;
                    while (y != y2)
                    {
                        first = horizontal13 ? static_cast<int>(xf + 0.5 * m13) : static_cast<int>(xf);
                        last = horizontal12 ? static_cast<int>(xl + 0.5 * m12) + 1 : static_cast<int>(xl) + 1;

                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - first) * step;

                        render_scanline<Shader, VertexData>(y,
                                                            first,
                                                            last,
                                                            shader.interpolate(v1, v3, w1f),
                                                            shader.interpolate(v1, v2, w1l),
                                                            interpolateF(ndc_v1[2], ndc_v3[2], w1f),
                                                            interpolateF(ndc_v1[2], ndc_v2[2], w1l),
                                                            w0,
                                                            step,
                                                            shader,
                                                            threadPool);
                        ++y;
                        if (y == height)
                        {
                            return;
                        }
                        w1f -= idy13;
                        w1l -= idy12;
                        xf += m13;
                        xl += m12;
                    }
                    if (y2 == y3)
                    {
                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - x3) * step;
                        render_scanline<Shader, VertexData>(y,
                                                            x3,
                                                            x2 + 1,
                                                            v3,
                                                            v2,
                                                            ndc_v3[2],
                                                            ndc_v2[2],
                                                            w0,
                                                            step,
                                                            shader,
                                                            threadPool);
                        return;
                    }
                    else
                    {
                        first = horizontal13 ? static_cast<int>(xf + 0.5 * m13) : static_cast<int>(xf);
                        step = 1.0f / (x2 + 1 - xf);
                        w0 = 1.0f + (xf - first) * step;

                        render_scanline<Shader, VertexData>(y,
                                                            first,
                                                            x2 + 1,
                                                            shader.interpolate(v1, v3, w1f),
                                                            v2,
                                                            interpolateF(ndc_v1[2], ndc_v3[2], w1f),
                                                            ndc_v2[2],
                                                            w0,
                                                            step,
                                                            shader,
                                                            threadPool);
                    }
                    ++y;
                    if (y == height)
                    {
                        return;
                    }
                    w1f -= idy13;
                    double w2l = (y3f - y) * idy23;
                    xl = m23 * y + q23;
                    xf += m13;

                    while (y != y3)
                    {
                        first = horizontal13 ? static_cast<int>(xf + 0.5 * m13) : static_cast<int>(xf);
                        last = horizontal23 ? static_cast<int>(xl + 0.5 * m23) + 1 : static_cast<int>(xl) + 1;

                        step = 1.0f / (xl - xf);
                        w0 = 1.0f + (xf - first) * step;

                        render_scanline<Shader, VertexData>(y, //-> 10
                                                            first,
                                                            last,
                                                            shader.interpolate(v1, v3, w1f),
                                                            shader.interpolate(v2, v3, w2l),
                                                            interpolateF(ndc_v1[2], ndc_v3[2], w1f),
                                                            interpolateF(ndc_v2[2], ndc_v3[2], w2l),
                                                            w0,
                                                            step,
                                                            shader,
                                                            threadPool);
                        ++y; // -> 11
                        if (y == height)
                        {
                            return;
                        }
                        w1f -= idy13;
                        w2l -= idy23;
                        xl += m23;
                        xf += m13;
                    }
                    step = 1.0f / (xl - xf);
                    w0 = 1.0f + (xf - x3) * step;
                    last = horizontal23 ? static_cast<int>(xl + 0.5 * m23) + 1 : static_cast<int>(xl) + 1;

                    render_scanline<Shader, VertexData>(y,
                                                        x3,
                                                        last,
                                                        v3,
                                                        shader.interpolate(v2, v3, w2l),
                                                        ndc_v3[2],
                                                        interpolateF(ndc_v2[2], ndc_v3[2], w2l),
                                                        w0,
                                                        step,
                                                        shader,
                                                        threadPool);
                }
            }
        }

        std::array<double, 16> projection_matrix;

    private:
        static inline double ndc2idxf(double ndc, int range)
        {
            return (ndc + 1.0f) * (range - 1) / 2.0f;
        }

        template <class VertexData>
        inline void project(const Vertex<VertexData> &v, std::array<double, 3> &ndc)
        {
            const double
                w = v.x * projection_matrix[4 * 3 + 0] + v.y * projection_matrix[4 * 3 + 1] + v.z * projection_matrix[4 * 3 + 2] + projection_matrix[4 * 3 + 3];
            ndc[0] = (v.x * projection_matrix[4 * 0 + 0] + v.y * projection_matrix[4 * 0 + 1] + v.z * projection_matrix[4 * 0 + 2] + projection_matrix[4 * 0 + 3]) / w;
            ndc[1] = (v.x * projection_matrix[4 * 1 + 0] + v.y * projection_matrix[4 * 1 + 1] + v.z * projection_matrix[4 * 1 + 2] + projection_matrix[4 * 1 + 3]) / w;
            ndc[2] = (v.x * projection_matrix[4 * 2 + 0] + v.y * projection_matrix[4 * 2 + 1] + v.z * projection_matrix[4 * 2 + 2] + projection_matrix[4 * 2 + 3]) / w;
        }

        static inline double interpolateF(double v1, double v2, double w)
        {
            return v1 * w + v2 * (1.0 - w);
        }

        /**
         * @brief The following method is the single-threaded scanline method.
         *
         * @tparam Shader
         * @tparam VertexData
         * @param pixel_y
         * @param xl
         * @param xr
         * @param vl
         * @param vr
         * @param ndczl
         * @param ndczr
         * @param w
         * @param step
         * @param shader
         */
        template <class Shader, class VertexData>
        void render_scanline(const int &pixel_y,
                             const int &xl,
                             const int &xr,
                             const Vertex<VertexData> &vl,
                             const Vertex<VertexData> &vr,
                             const double &ndczl,
                             const double &ndczr,
                             double w,
                             double step,
                             Shader &shader)
        {
            const double epsilon = 1.0e-8f;
            if (pixel_y < 0 || pixel_y >= height)
            {
                return;
            }
            int x = std::max(xl, 0);
            w += (xl - x) * step;
            for (; x != std::min(width, xr + 1); ++x)
            {
                const double ndcz = interpolateF(ndczl, ndczr, w);
                if ((z_buffer[pixel_y * width + x] + epsilon) < ndcz)
                {
                    continue;
                }
                Vertex<VertexData> p = shader.interpolate(vl, vr, w);
                shader.perspectiveCorrect(p);
                target[pixel_y * width + x] = shader.shade(p);
                z_buffer[pixel_y * width + x] = ndcz;
                w -= step;
            }
        }

        /**
         * @brief The following is the multi-threaded scanline method, which enforces access control to the zbuffer and targetT,
         * avoiding any buffer corruption.
         *
         * @tparam Shader
         * @tparam VertexData
         * @param pixel_y
         * @param xl
         * @param xr
         * @param vl
         * @param vr
         * @param ndczl
         * @param ndczr
         * @param w
         * @param step
         * @param shader
         * @param threadPool
         */
        template <class Shader, class VertexData>
        void render_scanline(const int &pixel_y,
                             const int &xl,
                             const int &xr,
                             const Vertex<VertexData> &vl,
                             const Vertex<VertexData> &vr,
                             const double &ndczl,
                             const double &ndczr,
                             double w,
                             double step,
                             Shader &shader,
                             thread_pool_manager<std::mutex *> &threadPool)
        {
            if (pixel_y < 0 || pixel_y >= height || xl > xr)
            {
                return;
            }
            int x = std::max(xl, 0);
            int index;
            double ndcz;
            Vertex<VertexData> p;
            w += (xl - x) * step;
            for (; x != std::min(width, xr + 1); ++x)
            {
                index = pixel_y * width + x;
                ndcz = interpolateF(ndczl, ndczr, w);
                {
                    std::lock_guard<std::mutex> lock(threadPool.protection[index]);
                    if ((z_buffer[index] + 1.0e-8f) < ndcz)
                    {
                        continue;
                    }
                    p = shader.interpolate(vl, vr, w);
                    shader.perspectiveCorrect(p);
                    target[index] = shader.shade(p);
                    z_buffer[index] = ndcz;
                }
                w -= step;
            }
        }

        int width;
        int height;
        ZBuffer_t *z_buffer;
        Target_t *target;
    };

    template <class Target_t>
    struct Object
    {
        Object() = default;

        virtual void render(Rasterizer<Target_t, double> &rasterizer,
                            thread_pool_manager<std::mutex *> &threadsPool,
                            bool isObjectLevelMultithreaded,
                            bool isPolygonLevelMultithreaded,
                            bool isRenderScanlineLevelMultithreaded) = 0;
    };

    template <class Target_t, class ShaderType, class VertexData>
    struct ObjectTemplate : Object<Target_t>
    {

        std::vector<Triangle3D<VertexData>> triangles;
        ShaderType shader;

        ObjectTemplate(
            ShaderType &&shader,
            std::vector<Triangle3D<VertexData>> &&triangles) : shader(std::move(shader)),
                                                               triangles(std::move(triangles)) {}

        /**
         * @brief The following method renders each triangle that makes up the object, allowing it in a single or multi threaded way.
         *
         * @param rasterizer
         * @param threadsPool
         */
        void render(Rasterizer<Target_t, double> &rasterizer,
                    thread_pool_manager<std::mutex *> &threadsPool,
                    bool isObjectLevelMultithreaded,
                    bool isPolygonLevelMultithreaded,
                    bool isRenderScanlineLevelMultithreaded)
        {
            if (isPolygonLevelMultithreaded)
            {
                if (isRenderScanlineLevelMultithreaded)
                {
                    // Polygon and Scanline
                    for (Triangle3D<VertexData> &triangle : triangles)
                    {
                        threadsPool.executeTask([&]
                                                { rasterizer.render_triangle(triangle.v0,
                                                                             triangle.v1,
                                                                             triangle.v2,
                                                                             shader,
                                                                             threadsPool,
                                                                             true); });
                    }
                }
                else
                {
                    // Polygon
                    for (Triangle3D<VertexData> &triangle : triangles)
                    {
                        threadsPool.executeTask([&]
                                                { rasterizer.render_triangle(triangle.v0,
                                                                             triangle.v1,
                                                                             triangle.v2,
                                                                             shader,
                                                                             threadsPool,
                                                                             false); });
                    }
                }
            }
            else
            {
                if (isRenderScanlineLevelMultithreaded)
                {
                    // Scanline
                    for (Triangle3D<VertexData> &triangle : triangles)
                    {
                        rasterizer.render_triangle(triangle.v0,
                                                   triangle.v1,
                                                   triangle.v2,
                                                   shader,
                                                   threadsPool,
                                                   true);
                    }
                }
                else
                {
                    // No Polygon and No ScanLine THREADS
                    if (isObjectLevelMultithreaded)
                    {
                        // Object
                        for (Triangle3D<VertexData> &triangle : triangles)
                        {
                            rasterizer.render_triangle(triangle.v0,
                                                       triangle.v1,
                                                       triangle.v2,
                                                       shader,
                                                       threadsPool,
                                                       false);
                        }
                    }
                    else
                    {
                        // NO LEVEL
                        for (Triangle3D<VertexData> &triangle : triangles)
                        {
                            rasterizer.render_triangle(triangle.v0,
                                                       triangle.v1,
                                                       triangle.v2,
                                                       shader);
                        }
                    }
                }
            }
        }
    };

    template <class Target_t>
    struct Scene
    {
        /**
         * @brief The following vector contains all the objects that make up the scene.
         *
         */
        std::vector<Object<Target_t> *> objects;

        /**
         * @brief The following constructor takes all the objects needed to compose the scene.
         *
         * @param objects
         */
        Scene(std::vector<Object<Target_t> *> &&objects) : objects(std::move(objects)) {}

        /**
         * @brief The following method renders each scene's object, allowing for single or multi threaded execution.
         *
         * @param rasterizer
         * @param threadsPool
         * @param isObjectLevelMultithreaded
         * @param isPolygonLevelMultithreaded
         * @param isRenderScanlineLevelMultithreaded
         */
        void render(Rasterizer<Target_t, double> &rasterizer,
                    thread_pool_manager<std::mutex *> &threadsPool,
                    bool isObjectLevelMultithreaded,
                    bool isPolygonLevelMultithreaded,
                    bool isRenderScanlineLevelMultithreaded)
        {
            if (isObjectLevelMultithreaded)
            {
                for (Object<Target_t> *&object : objects)
                {
                    threadsPool.executeTask([&]
                                            { object->render(
                                                  rasterizer,
                                                  threadsPool,
                                                  isObjectLevelMultithreaded,
                                                  isPolygonLevelMultithreaded,
                                                  isRenderScanlineLevelMultithreaded); });
                }
            }
            else
            {
                for (Object<Target_t> *&object : objects)
                {
                    object->render(
                        rasterizer,
                        threadsPool,
                        isObjectLevelMultithreaded,
                        isPolygonLevelMultithreaded,
                        isRenderScanlineLevelMultithreaded);
                }
            }
        }
    };

}; // namespace pipeline3D

#endif // RASTERIZATION_HPP

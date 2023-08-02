#include <iostream>
#include <chrono>
#include "rasterization.hpp"

using namespace pipeline3D;

struct TextureData
{
	double u;
	double v;

	TextureData() = default;

	TextureData(double u, double v) : u(u), v(v) {}
};

struct TextureShader : public VertexShader<Pixel, TextureData, TextureShader>
{
	BitmapImage texture;

	TextureShader(BitmapImage &&texture) : texture(std::forward<BitmapImage>(texture)) {}

	/**
	 * @brief The following functions defines the interpolation needed for vertices whose extra data is TextureData.
	 *
	 * @param v1
	 * @param v2
	 * @param w
	 * @return Vertex<TextureData>
	 */
	Vertex<TextureData> interpolate(const Vertex<TextureData> &v1,
									const Vertex<TextureData> &v2,
									double w)
	{
		const double w2 = (1.0f - w);
		return Vertex<TextureData>(
			(w * v1.x + w2 * v2.x), // x,
			(w * v1.y + w2 * v2.y), // y,
			(w * v1.z + w2 * v2.z), // z,
			TextureData({
				v1.vertexData.u * w + v2.vertexData.u * w2,
				v1.vertexData.v * w + v2.vertexData.v * w2,
			})); // VertexData
	}

	/**
	 * @brief The following function defines the perspective correction needed for vertices whose extra data is TextureData.
	 *
	 * @param v
	 */
	void perspectiveCorrect(Vertex<TextureData> &v)
	{
		v.z = 1.0f / v.z;
		v.x *= v.z;
		v.y *= v.z;
		v.vertexData.u *= v.z;
		v.vertexData.v *= v.z;
	}

	Pixel shade(Vertex<TextureData> v)
	{
		Pixel got = texture.pixel_for_uv(v.vertexData.u, v.vertexData.v);
		return {got.r, got.g, got.b};
	}
};

struct PixelShader : public VertexShader<Pixel, Pixel, PixelShader>
{
	/**
	 * @brief The following functions defines the interpolation needed for vertices whose extra data
	 * is Pixel, meaning that it will be used to interpolate colored vertices.
	 *
	 */
	Vertex<Pixel> interpolate(const Vertex<Pixel> &v1,
							  const Vertex<Pixel> &v2,
							  double w)
	{
		const double w2 = (1.0f - w);
		return Vertex<Pixel>(
			(w * v1.x + w2 * v2.x), // x,
			(w * v1.y + w2 * v2.y), // y,
			(w * v1.z + w2 * v2.z), // z,
			Pixel({static_cast<uint8_t>(v1.vertexData.r * w + v2.vertexData.r * w2),
				   static_cast<uint8_t>(v1.vertexData.g * w + v2.vertexData.g * w2),
				   static_cast<uint8_t>(v1.vertexData.b * w + v2.vertexData.b * w2)})); // VertexData
	}

	/**
	 * @brief The following function defines the perspective correction needed for vertices whose extra data is Pixel, meaning
	 * that it will be used to interpolate colored vertices.
	 *
	 * @param v
	 */
	void perspectiveCorrect(Vertex<Pixel> &v)
	{
		v.z = 1.0f / v.z;
		v.x *= v.z;
		v.y *= v.z;
	}

	Pixel shade(Vertex<Pixel> v)
	{
		return v.vertexData;
	}
};

int main()
{
	const bool benchmark = false;
	const bool isObjectLevelMultithreaded = false;
	const bool isPolygonLevelMultithreaded = false;
	const bool isRenderScanlineLevelMultithreaded = true;
	constexpr bool isMultithreaded = isObjectLevelMultithreaded || isPolygonLevelMultithreaded || isRenderScanlineLevelMultithreaded;
	const bool writeOut = true;
	// Working formats: 1000, 4000, 5000, 6000, 10000
	const int w = 8000;
	const int h = 8000;
	const int fieldSize = w * h;
	constexpr int numberOfThreads = 25;
	Rasterizer<Pixel, double> rasterizer;

	rasterizer.set_perspective_projection(-1, 1, -1, 1, 1, 2);
	std::vector<Pixel> target(fieldSize, Pixel({0, 0, 0}));
	std::vector<double> zBuffer(fieldSize, 1.0e8f);
	rasterizer.set_target(w, h, &target[0], &zBuffer[0]);

	constexpr Pixel red({255, 0, 0});
	constexpr Pixel green({0, 255, 0});
	constexpr Pixel blue({0, 0, 255});
	constexpr Pixel yellow({255, 255, 0});
	constexpr float slope = -0.2f;

	// gradient
	const Vertex<Pixel>
		v1 = {1.2, -0.5, 2.5f + slope * (1 - 1), blue},
		v2 = {1, 1, 1.5f + slope * (1 + 1), yellow},
		v3 = {-0.2, 1, 2.5f - +slope * (-1 + 1), red},
		v4 = {0.2, 0.2, 1.5f + slope * (-1 - 1), green};

	// neon
	const Vertex<TextureData>
		d1 = {-1, -1.2, 1.5f + slope * (1 - 1), {0, 0}},
		d2 = {-1, 0.2, 1.5f + slope * (1 + 1), {0, 1}},
		d3 = {0.1, -0., 1.5f - +slope * (-1 + 1), {1, 1}},
		d4 = {-0., -1.3, 1.5f + slope * (-1 - 1), {1, 0}};

	// bear
	const Vertex<TextureData>
		p1 = {0.2, -1, 1.5f + slope * (1 - 1), {0, 0}},
		p2 = {-0., 0.3, 1.5f + slope * (1 + 1), {0, 1}},
		p3 = {1.3, 0.2, 1.5f - +slope * (-1 + 1), {1, 1}},
		p4 = {1.4, -1, 1.5f + slope * (-1 - 1), {1, 0}};

	// heinstein
	const Vertex<TextureData>
		l1 = {-1, 0., 1.5f + slope * (1 - 1), {0, 0}},
		l2 = {-1, 1, 1.5f + slope * (1 + 1), {0, 1}},
		l3 = {0.15, 1, 1.5f - +slope * (-1 + 1), {1, 1}},
		l4 = {0., 0.2, 1.5f + slope * (-1 - 1), {1, 0}};

	// bear 2 (overlapped by heinstein)
	const Vertex<TextureData>
		m1 = {-1, 0., 2.5f + slope * (1 - 1), {0, 0}},
		m2 = {-1, 1, 2.5f + slope * (1 + 1), {0, 1}},
		m3 = {0.15, 1, 1.5f - +slope * (-1 + 1), {1, 1}},
		m4 = {0., 0.2, 1.5f + slope * (-1 - 1), {1, 0}};

	ObjectTemplate<Pixel, TextureShader, TextureData> heinsteinMuralesObject(
		//!--------------------EACH PATHNAME needs to be overwritten with the path where the images are stored.-----------------------------------
		TextureShader(BitmapImage("............................./heinsteinMurales.bmp")),
		{
			{l1, l2, l3},
			{l4, l1, l3},
		});

	ObjectTemplate<Pixel, TextureShader, TextureData> neonWorldObject(
		//!--------------------EACH PATHNAME needs to be overwritten with the path where the images are stored.-----------------------------------
		TextureShader(BitmapImage("............................./neonWorld4K.bmp")),
		{
			{d1, d2, d3},
			{d4, d1, d3},
		});

	ObjectTemplate<Pixel, PixelShader, Pixel> gradientObject(
		PixelShader(),
		{
			{v1, v2, v3},
			{v4, v1, v3},
		});

	ObjectTemplate<Pixel, TextureShader, TextureData> bearObject(
		//!--------------------EACH PATHNAME needs to be overwritten with the path where the images are stored.-----------------------------------
		TextureShader(BitmapImage("............................./bear8k.bmp")),
		{
			{p1, p2, p3},
			{p4, p1, p3},
		});

	ObjectTemplate<Pixel, TextureShader, TextureData> bearObject2(
		//!--------------------EACH PATHNAME needs to be overwritten with the path where the images are stored.-----------------------------------
		TextureShader(BitmapImage("............................./bear8k.bmp")),
		{
			{m1, m2, m3},
			{m4, m1, m3},
		});

	Scene<Pixel> scene({&gradientObject, &bearObject, &neonWorldObject, &bearObject2, &heinsteinMuralesObject});

	if (benchmark)
	{
		const uint8_t testCycles = 30;
		std::cout << "starting benchmark for " << w << "x" << h << " target" << std::endl;
		std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
		// 25 Threads perform better with a Windows10 + Intel i7_10875H mobile (8 core/ 16 Threads)
		thread_pool_manager<std::mutex *> pool(numberOfThreads);
		if (isMultithreaded)
		{
			std::cout << "Number of threads created: " << pool.getPoolSize() << std::endl;
			pool.protection = new std::mutex[fieldSize];
			start_time = std::chrono::high_resolution_clock::now();
			for (uint8_t i = 0; i < testCycles; i++)
			{
				scene.render(rasterizer, pool, isObjectLevelMultithreaded, isPolygonLevelMultithreaded, isRenderScanlineLevelMultithreaded);
			}
			pool.shutdown();
			// protection method deallocation.
			delete[] pool.protection;
			std::chrono::time_point<std::chrono::high_resolution_clock> end_time = std::chrono::high_resolution_clock::now();
			float elapsed_time = std::chrono::duration<float>(end_time - start_time).count();
			std::cout << "elapsed time: " << elapsed_time << std::endl;
			std::cout << "Done Rendering" << std::endl;
			// final print.
			if (writeOut)
			{
				//!--------------------This pathnames must be overwritten with the desired output path.-----------------------------------
				stbi_write_bmp("............................./out_scene.bmp", w, h, 3, &target[0]);
			}
		}
		else
		{
			pool.shutdown();
			start_time = std::chrono::high_resolution_clock::now();
			for (uint8_t i = 0; i < testCycles; i++)
			{
				scene.render(rasterizer, pool, isObjectLevelMultithreaded, isPolygonLevelMultithreaded, isRenderScanlineLevelMultithreaded);
			}
			std::chrono::time_point<std::chrono::high_resolution_clock> end_time = std::chrono::high_resolution_clock::now();
			float elapsed_time = std::chrono::duration<float>(end_time - start_time).count();
			std::cout << "elapsed time: " << elapsed_time << std::endl;
			std::cout << "Done Rendering" << std::endl;
			// final print.
			if (writeOut)
			{
				//!--------------------This pathnames must be overwritten with the desired output path.-----------------------------------
				stbi_write_bmp("............................./out_scene.bmp", w, h, 3, &target[0]);
			}
		}
	}
	else
	{
		// 25 Threads perform better with a Windows10 + Intel i7_10875H mobile (8 core/ 16 Threads)
		thread_pool_manager<std::mutex *> pool(numberOfThreads);
		std::cout << "starting benchmark for " << w << "x" << h << " target" << std::endl;
		std::cout << "Write on file is: " << writeOut << std::endl;
		std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
		if (isMultithreaded)
		{
			std::cout << "Number of threads created: " << pool.getPoolSize() << std::endl;
			pool.protection = new std::mutex[fieldSize];
			start_time = std::chrono::high_resolution_clock::now();
			scene.render(rasterizer, pool, isObjectLevelMultithreaded, isPolygonLevelMultithreaded, isRenderScanlineLevelMultithreaded);
			pool.shutdown();
			// protection method deallocation.
			delete[] pool.protection;
			std::chrono::time_point<std::chrono::high_resolution_clock> end_time = std::chrono::high_resolution_clock::now();
			float elapsed_time = std::chrono::duration<float>(end_time - start_time).count();
			std::cout << "elapsed time: " << elapsed_time << std::endl;
			std::cout << "Done Rendering" << std::endl;
			// final print.
			if (writeOut)
			{
				//!--------------------This pathnames must be overwritten with the desired output path.-----------------------------------
				stbi_write_bmp("............................./out_scene.bmp", w, h, 3, &target[0]);
			}
		}
		else
		{
			pool.shutdown();
			start_time = std::chrono::high_resolution_clock::now();
			scene.render(rasterizer, pool, isObjectLevelMultithreaded, isPolygonLevelMultithreaded, isRenderScanlineLevelMultithreaded);
			auto end_time = std::chrono::high_resolution_clock::now();
			float elapsed_time = std::chrono::duration<float>(end_time - start_time).count();
			std::cout << "elapsed time: " << elapsed_time << std::endl;
			std::cout << "Done Rendering" << std::endl;
			// final print.
			if (writeOut)
			{
				//!--------------------This pathnames must be overwritten with the desired output path.-----------------------------------
				stbi_write_bmp("............................./out_scene.bmp", w, h, 3, &target[0]);
			}
		}
	}

	return 0;
}

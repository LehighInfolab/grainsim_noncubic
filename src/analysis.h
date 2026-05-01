//analysis.h
#include <iostream>
#include <fstream>
#include <unordered_map>

#include "types.h"
#include "lattice.h"


#include "json.hpp"
using json = nlohmann::json;


class lattice_analyzer_t
{
private:

	lattice_t *curr_cube;

	spin_t max_grains;
	spin_t calculate_max_grains()
	{
		spin_t max_so_far = 0;
		for (coord_t z = 0; z < curr_cube->side_length_z; ++z)
			for (coord_t y = 0; y < curr_cube->side_length_y; ++y)
				for (coord_t x = 0; x < curr_cube->side_length_x; ++x)
				{
					spin_t curr_spin = curr_cube->voxel_at(x, y, z)->spin;
					if (curr_spin > max_so_far) max_so_far = curr_spin;
				}
		return max_so_far;
	}

#pragma pack(push, 1)
	struct boundary_info_t
	{
		int sm_to_lg_outies = 0,
			sm_to_lg_delta = 0,
			lg_to_sm_outies = 0,
			lg_to_sm_delta = 0,
			surface_area = 0;
	};
#pragma pack(pop)
	std::unordered_map<spin_t, std::unordered_map<spin_t, boundary_info_t> > sparse_info_matrix;
	std::unordered_map<spin_t, size_t> vol_map;

	void incr_sparse_outies(spin_t a, spin_t b)
	{
		if (a > b)
		{
			++sparse_info_matrix[b][a].lg_to_sm_outies;
		}
		else if (a < b)
		{
			++sparse_info_matrix[a][b].sm_to_lg_outies;
		}
	}
	void incr_sparse_delta(spin_t a, spin_t b)
	{
		if (a > b)
		{
			++sparse_info_matrix[b][a].lg_to_sm_delta;
		}
		else if (a < b)
		{
			++sparse_info_matrix[a][b].sm_to_lg_delta;
		}
	}
	void incr_sparse_sa(spin_t a, spin_t b)
	{
		++sparse_info_matrix[a > b ? b : a][a > b ? a : b].surface_area;
	}

	void check_edge(
		coord_t rx, coord_t ry, coord_t rz,
		coord_t x1, coord_t y1, coord_t z1,
		coord_t x2, coord_t y2, coord_t z2,
		coord_t x3, coord_t y3, coord_t z3,
		coord_t x4, coord_t y4, coord_t z4)
	{
		spin_t id1, id2, id3, id4;

		id1 = curr_cube->voxel_at(rx + x1, ry + y1, rz + z1)->spin;
		id2 = curr_cube->voxel_at(rx + x2, ry + y2, rz + z2)->spin;
		id3 = curr_cube->voxel_at(rx + x3, ry + y3, rz + z3)->spin;
		id4 = curr_cube->voxel_at(rx + x4, ry + y4, rz + z4)->spin;

		if (id1 != id2 && id2 == id3 && id2 == id4)
		{
			incr_sparse_outies(id1, id2);
		}
		else if (id2 != id1 && id1 == id3 && id1 == id4)
		{
			incr_sparse_outies(id2, id1);
		}
		else if (id3 != id1 && id1 == id2 && id1 == id4)
		{
			incr_sparse_outies(id3, id1);
		}
		else if (id4 != id1 && id1 == id2 && id1 == id3)
		{
			incr_sparse_outies(id4, id1);
		}
	}

	size_t matrix_dim = 0;
	void generate_matrices()
	{
		sparse_info_matrix.clear();
		vol_map.clear();

		for(coord_t z = 0; z < curr_cube->side_length_z; ++z)
			for (coord_t y = 0; y < curr_cube->side_length_y; ++y)
				for (coord_t x = 0; x < curr_cube->side_length_x; ++x)
				{
					spin_t
						curr_id = curr_cube->voxel_at(x, y, z)->spin,
						fwd_id = curr_cube->voxel_at(x, y, z + 1)->spin,
						right_id = curr_cube->voxel_at(x + 1, y, z)->spin,
						up_id = curr_cube->voxel_at(x, y + 1, z)->spin;

					if (curr_id != right_id)
					{
						incr_sparse_sa(curr_id, right_id);
					}
					if (curr_id != fwd_id)
					{
						incr_sparse_sa(curr_id, fwd_id);
					}
					if (curr_id != up_id)
					{
						incr_sparse_sa(curr_id, up_id);
					}

					// back bottom
					check_edge(
						x, y, z,
						0, 0, -1, 
						0, 0, 0,
						0, -1, 0,
						0, -1, -1);
					// back left
					check_edge(
						x, y, z,
						-1, 0, 0,
						0, 0, 0, 
						0, 0, -1, 
						-1, 0, -1);
					// top left
					check_edge(
						x, y, z,
						-1, 1, 0, 
						0, 1, 0,
						0, 0, 0,
						-1, 0, 0);

					++vol_map[curr_id];
				}

	}
	// calculate all centroid
	struct centroid_t {
		double x= 0, y=0, z=0;
	};

	// compute all centroids
	std::unordered_map<spin_t, centroid_t> compute_centroids() const 
	{
		std::unordered_map<spin_t, centroid_t> sums;
		sums.reserve(vol_map.size());

		// sum each grain's information
		for (coord_t z=0; z<curr_cube->side_length_z; ++z) {
			for (coord_t y = 0; y < curr_cube->side_length_y; ++y) {
				for (coord_t x = 0; x < curr_cube->side_length_x; ++x) {
					spin_t id = curr_cube->voxel_at(x, y, z)->spin;
					sums[id].x += static_cast<double>(x);
					sums[id].y += static_cast<double>(y);
					sums[id].z += static_cast<double>(z);
				}
			}
		}

		std::unordered_map<spin_t, centroid_t> centroids;
		centroids.reserve(vol_map.size());

		for (const auto& [id, vol] : vol_map) 
		{
			if (vol==0) continue;
			double v = static_cast<double>(vol);
			auto it = sums.find(id);
			if (it != sums.end()) {
				centroids[id] = {it->second.x / v, it->second.y / v, it->second.z / v};
			}
		}
		return centroids;
	}

	std::unordered_map<spin_t, double> compute_zcenters_voxel_index() const
    {
        std::unordered_map<spin_t, double> z_sum;
        z_sum.reserve(vol_map.size());

        for (coord_t z = 0; z < curr_cube->side_length_z; ++z)
            for (coord_t y = 0; y < curr_cube->side_length_y; ++y)
                for (coord_t x = 0; x < curr_cube->side_length_x; ++x)
                {
                    spin_t id = curr_cube->voxel_at(x, y, z)->spin;
                    z_sum[id] += static_cast<double>(z);
                }

        std::unordered_map<spin_t, double> zcenter;
        zcenter.reserve(vol_map.size());

        for (const auto& [id, vol] : vol_map)
        {
            if (vol == 0) continue;
            auto it = z_sum.find(id);
            double sum = (it == z_sum.end()) ? 0.0 : it->second;
            zcenter[id] = sum / static_cast<double>(vol);
        }

        return zcenter;
    }

public:
	size_t get_grain_count() const {
		return vol_map.size();
	}

	void load_lattice(lattice_t *cube)
	{
		curr_cube = cube;
		max_grains = calculate_max_grains();
		matrix_dim = max_grains + 1;

		generate_matrices();
	}

	double get_curvature(spin_t a, spin_t b) // DOES NOT VERIFY THAT BOUNDARY EXISTS!!!!!
	{
		if (a > b)
		{
			return (3.141592653589793 / 4.0) * (sparse_info_matrix[b][a].lg_to_sm_outies - sparse_info_matrix[b][a].sm_to_lg_outies);
		}
		else if (a < b)
		{
			return (3.141592653589793 / 4.0) * (sparse_info_matrix[a][b].sm_to_lg_outies - sparse_info_matrix[a][b].lg_to_sm_outies);
		}
		return 0;
	}

	void save_analysis_to_file(const char *fname)
	{
		std::ofstream afile(fname);

		std::cout << "Creating analysis file " << fname << std::endl;

		// z center
		const auto zcenter = compute_zcenters_voxel_index();
		// Volumes
		afile << "VOLUMES\n";
		// for (auto vol_iter = vol_map.begin(); vol_iter != vol_map.end(); ++vol_iter)
		// {
		// 	afile << vol_iter->first << ' ' << vol_iter->second << '\n';
		// }

		for (const auto& [id, vol] : vol_map)
        {
            double cz = 0.0;
            auto it = zcenter.find(id);
            if (it != zcenter.end()) cz = it->second;

            afile << id << ' ' << vol << ' ' << cz << '\n';
        }
		
		// Curvatures
		afile << "CURVATURES\n";
		for (auto sm_iter = sparse_info_matrix.begin(); sm_iter != sparse_info_matrix.end(); ++sm_iter)
		{
			for (auto lg_iter = sm_iter->second.begin(); lg_iter != sm_iter->second.end(); ++lg_iter)
			{
				if (lg_iter->second.surface_area == 0) continue;

				afile << sm_iter->first << ' ' << lg_iter->first << ' ' << get_curvature(sm_iter->first, lg_iter->first) << '\n';
				afile << lg_iter->first << ' ' << sm_iter->first << ' ' << get_curvature(lg_iter->first, sm_iter->first) << '\n';
			}
		}

		// Curvatures
		afile << "SURFACE_AREAS\n";
		for (auto sm_iter = sparse_info_matrix.begin(); sm_iter != sparse_info_matrix.end(); ++sm_iter)
		{
			for (auto lg_iter = sm_iter->second.begin(); lg_iter != sm_iter->second.end(); ++lg_iter)
			{
				if (lg_iter->second.surface_area == 0) continue;

				afile << sm_iter->first << ' ' << lg_iter->first << ' ' << lg_iter->second.surface_area << '\n';
				afile << lg_iter->first << ' ' << sm_iter->first << ' ' << lg_iter->second.surface_area << '\n';
			}
		}

		// Velocities
		afile << "VELOCITIES\n";
		for (auto sm_iter = curr_cube->boundary_tracker.velocity_tracker.begin(); sm_iter != curr_cube->boundary_tracker.velocity_tracker.end(); ++sm_iter)
		{
			for (auto lg_iter = sm_iter->second.begin(); lg_iter != sm_iter->second.end(); ++lg_iter)
			{
				std::pair<int, int> delta = lg_iter->second;
				
				afile << sm_iter->first << ' ' << lg_iter->first << ' ' << (delta.first - delta.second) << '\n';
				afile << lg_iter->first << ' ' << sm_iter->first << ' ' << (delta.second - delta.first) << '\n';
			}
		}

		afile << "ADJACENT_BOUNDARIES\n";
		for (auto sm_iter = curr_cube->boundary_tracker.boundary_map.begin(); sm_iter != curr_cube->boundary_tracker.boundary_map.end(); ++sm_iter)
		{
			for (auto lg_iter = sm_iter->second.begin(); lg_iter != sm_iter->second.end(); ++lg_iter)
			{
				boundary_t *boundary = lg_iter->second;

				if (boundary->area() == 0) continue;

				afile << boundary->a_spin << '/' << boundary->b_spin;
				for (auto junc_iter = boundary->junctions.begin(); junc_iter != boundary->junctions.end(); ++junc_iter)
				{
					afile << ' ' << junc_iter->first->a_spin << '/' << junc_iter->first->b_spin;
				}

				afile << '\n';
			}
		}

		curr_cube->boundary_tracker.reset_flip_tracker();

		afile.close();
	}

	void save_analysis_to_json(const char *fname, double timestep) {
		std::cout << "Creating Json analysis file " << fname << std::endl;

		const auto centroids = compute_centroids();

		struct boundary_data_t {
			spin_t neighbor;
			double curvature;
			int surface_area;
			int velocity;
			double mobility;
			bool transformed;
		};
		std::unordered_map<spin_t, std::vector<boundary_data_t>> grain_boundaries;

		for (const auto& [a, inner] : sparse_info_matrix) {
			for (const auto& [b, info] : inner) {
				if (info.surface_area == 0) continue;

				double curv_a = get_curvature(a, b);
				double curv_b = get_curvature(b, a);

				// get Mobility and transformed info 
				bool trans = curr_cube->boundary_tracker.is_transformed(a, b);
				double mobility = trans ? curr_cube->transitioned_mobility : curr_cube->default_mobility;

				// Velocity
				int vel_ab = 0, vel_ba = 0;
				auto v_it = curr_cube->boundary_tracker.velocity_tracker.find(a);
				if (v_it != curr_cube->boundary_tracker.velocity_tracker.end()) {
					auto v_it2 = v_it->second.find(b);
					if (v_it2 != v_it->second.end()) {
						vel_ab = v_it2->second.first - v_it2->second.second;
						vel_ba = v_it2->second.second - v_it2->second.first;
					}
				}

				grain_boundaries[a].push_back({b, curv_a, info.surface_area, vel_ab, mobility, trans});
				grain_boundaries[b].push_back({a, curv_b, info.surface_area, vel_ba, mobility, trans});
			}
		}


		// construct json output
		json output;
		output["metadata"] = {
			{"timestamp_mcs", timestep},
			{"total_grains", vol_map.size()},
			{"lattice_dimensions", {curr_cube->side_length_x, curr_cube->side_length_y, curr_cube->side_length_z}}
		};

		json grains = json::array();

		for (const auto& [id, vol] : vol_map) {
			centroid_t c = {0, 0, 0};
			auto c_it = centroids.find(id);
			if (c_it != centroids.end()) c = c_it->second;

			json boundaries = json::array();
			json neighbors = json::array();

			auto b_it = grain_boundaries.find(id);
			if (b_it != grain_boundaries.end()) {
				for (const auto& bd : b_it->second) {
					neighbors.push_back(bd.neighbor);
					boundaries.push_back({
						{"neighbor_id",   bd.neighbor},
						{"curvature",     bd.curvature},
						{"surface_area",  bd.surface_area},
						{"velocity",      bd.velocity},
						{"mobility",      bd.mobility},
						{"transformed",   bd.transformed}
					});
				}
			}

			// v = (4 * pi r^3 )/3
			double eq_radius = std::cbrt((3.0 * vol) / (4.0 * 3.14));

			grains.push_back({
				{"id",                (int)id},
				{"volume",            (int)vol},
				{"equivalent_radius", eq_radius},
				{"centroid",          {{"x", c.x}, {"y", c.y}, {"z", c.z}}},
				{"neighbors",         neighbors},
				{"boundaries",        boundaries}
			});
		}

		output["grains"] = grains;

		std::ofstream afile(fname);
		afile << output.dump(2);
		afile.close();

		curr_cube->boundary_tracker.reset_flip_tracker();
	}

};
// zip_utils.h
// Utility functions for bundling simulation output files into zip archives.
// Uses per-checkpoint zipping with a final merge to balance disk usage and simplicity.
#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <filesystem>
#include <zip.h>

// Bundle a VTK file and an optional JSON file into a single zip archive,
// then delete the originals. The zip file name is derived from the VTK path
// by replacing .vtk with .zip.
inline void bundle_and_cleanup(const std::string& vtk_path, const std::string& json_path) {
    std::string zip_path = vtk_path.substr(0, vtk_path.size() - 4) + ".zip";
    int zerr = 0;
    zip_t* archive = zip_open(zip_path.c_str(), ZIP_CREATE | ZIP_TRUNCATE, &zerr);
    if (!archive) {
        std::cerr << "Failed to create zip: " << zip_path << std::endl;
        return;
    }

    // Add VTK file with max compression
    zip_source_t* src = zip_source_file(archive, vtk_path.c_str(), 0, -1);
    if (src) {
        auto idx = zip_file_add(archive, std::filesystem::path(vtk_path).filename().c_str(), src, ZIP_FL_OVERWRITE);
        zip_set_file_compression(archive, idx, ZIP_CM_DEFLATE, 9);
    }

    // Add JSON file if provided
    if (!json_path.empty()) {
        zip_source_t* jsrc = zip_source_file(archive, json_path.c_str(), 0, -1);
        if (jsrc) {
            auto idx = zip_file_add(archive, std::filesystem::path(json_path).filename().c_str(), jsrc, ZIP_FL_OVERWRITE);
            zip_set_file_compression(archive, idx, ZIP_CM_DEFLATE, 9);
        }
    }

    // zip_close is where the actual compression and writing happens.
    // Only delete originals if the zip was written successfully.
    if (zip_close(archive) == 0) {
        std::filesystem::remove(vtk_path);
        if (!json_path.empty()) std::filesystem::remove(json_path);
    } else {
        std::cerr << "zip_close failed for: " << zip_path << std::endl;
    }
}

// Merge multiple small zip files into one final zip archive,
// then delete the small zips. Uses ZIP_FL_COMPRESSED to copy
// already-compressed data directly without re-compressing,
// so this runs at roughly file-copy speed.
inline void merge_zips(const std::vector<std::string>& small_zips, const std::string& final_path) {
    int zerr = 0;
    zip_t* final_archive = zip_open(final_path.c_str(), ZIP_CREATE | ZIP_TRUNCATE, &zerr);
    if (!final_archive) {
        std::cerr << "Failed to create final zip: " << final_path << std::endl;
        return;
    }

    // Source zips must stay open until zip_close(final_archive),
    // because ZIP_FL_COMPRESSED stores references, not copies.
    std::vector<zip_t*> sources;

    for (const auto& zf : small_zips) {
        int e = 0;
        zip_t* sz = zip_open(zf.c_str(), ZIP_RDONLY, &e);
        if (!sz) continue;
        sources.push_back(sz);

        for (zip_int64_t i = 0; i < zip_get_num_entries(sz, 0); ++i) {
            zip_source_t* src = zip_source_zip_file(final_archive, sz, i, ZIP_FL_COMPRESSED, 0, -1, nullptr);
            if (src) {
                auto idx = zip_file_add(final_archive, zip_get_name(sz, i, 0), src, ZIP_FL_OVERWRITE);
                zip_set_file_compression(final_archive, idx, ZIP_CM_DEFLATE, 0);
            }
        }
    }

    // Actual writing happens here — reads compressed data from source zips
    bool ok = (zip_close(final_archive) == 0);

    // Now safe to close sources
    for (auto* s : sources) zip_close(s);

    if (ok) {
        for (const auto& zf : small_zips) std::filesystem::remove(zf);
        std::cout << "Created final zip: " << final_path << std::endl;
    } else {
        std::cerr << "merge zip_close failed" << std::endl;
    }
}
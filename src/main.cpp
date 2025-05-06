#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/repair.h> 
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/centroid.h>
#include <CGAL/Circle_3.h>
#include <CGAL/Plane_3.h>

#include <boost/program_options.hpp>
#include <filesystem>
#include <iostream>
#include <vector>
#include <string>
#include <cmath> 

// Define Kernel and Mesh types
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Point_3<Kernel> Point_3;
typedef CGAL::Vector_3<Kernel> Vector_3;
typedef CGAL::Iso_cuboid_3<Kernel> Iso_cuboid_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef CGAL::Aff_transformation_3<Kernel> Aff_transformation_3;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace CGALIO = CGAL::Polygon_mesh_processing::IO;
namespace po = boost::program_options;
namespace fs = std::filesystem;


Mesh make_triangulated_cube(const Point_3& center, const Vector_3& dimensions) {
    Mesh cube_mesh;
    FT hx = dimensions.x() / FT(2);
    FT hy = dimensions.y() / FT(2);
    FT hz = dimensions.z() / FT(2);

    std::vector<Point_3> points = {
        center + Vector_3(-hx, -hy, -hz), // 0: left bottom back
        center + Vector_3( hx, -hy, -hz), // 1: right bottom back
        center + Vector_3( hx,  hy, -hz), // 2: right top back
        center + Vector_3(-hx,  hy, -hz), // 3: left top back
        center + Vector_3(-hx, -hy,  hz), // 4: left bottom front
        center + Vector_3( hx, -hy,  hz), // 5: right bottom front
        center + Vector_3( hx,  hy,  hz), // 6: right top front
        center + Vector_3(-hx,  hy,  hz)  // 7: left top front
    };

    std::vector<typename Mesh::Vertex_index> vertices;
    for(const auto& p : points) {
        vertices.push_back(cube_mesh.add_vertex(p));
    }

    // Add faces (triangulated)
    // Bottom face (-Z)
    cube_mesh.add_face(vertices[0], vertices[1], vertices[2]);
    cube_mesh.add_face(vertices[0], vertices[2], vertices[3]);
    // Top face (+Z)
    cube_mesh.add_face(vertices[4], vertices[6], vertices[5]);
    cube_mesh.add_face(vertices[4], vertices[7], vertices[6]);
    // Back face (-Y)
    cube_mesh.add_face(vertices[0], vertices[4], vertices[5]);
    cube_mesh.add_face(vertices[0], vertices[5], vertices[1]);
    // Front face (+Y)
    cube_mesh.add_face(vertices[3], vertices[6], vertices[7]);
    cube_mesh.add_face(vertices[3], vertices[2], vertices[6]);
    // Left face (-X)
    cube_mesh.add_face(vertices[0], vertices[7], vertices[4]);
    cube_mesh.add_face(vertices[0], vertices[3], vertices[7]);
    // Right face (+X)
    cube_mesh.add_face(vertices[1], vertices[5], vertices[6]);
    cube_mesh.add_face(vertices[1], vertices[6], vertices[2]);

    return cube_mesh;
}


// Cylinder Generation (Using Circle_3 + Internal Repair) 
// NB: This was unusually tricky to get right, with many edge cases.
// Perhaps due to the num_segments being too high or low, or the cylinder being too small.
// The CGAL library had some issues with degenerate faces and non-manifold edges.
// Hence, a known working Cylinder STL is used for pegs/holes.
// See: ../assets/peg.stl
Mesh make_cylinder(const Point_3& base_center, const Point_3& top_center, const FT radius, const int num_segments = 32) {
    Mesh cylinder;

    Vector_3 main_direction = top_center - base_center;
    FT V_sq_len = main_direction.squared_length();

    if (V_sq_len < FT(1e-18)) {
        std::cerr << "Error in make_cylinder: Base and top centers are coincident." << std::endl;
        return Mesh();
    }

    // Normalize the main direction vector
    double V_len = std::sqrt(CGAL::to_double(V_sq_len));
    Vector_3 normal = main_direction / FT(V_len);

    // Define the base and top circles
    CGAL::Circle_3<Kernel> base_circle(base_center, radius, normal);
    CGAL::Circle_3<Kernel> top_circle(top_center, radius, normal);

    // Get the supporting planes and their basis vectors
    CGAL::Plane_3<Kernel> base_plane = base_circle.supporting_plane();
    CGAL::Plane_3<Kernel> top_plane = top_circle.supporting_plane();
    Vector_3 base_u = base_plane.base1(); Vector_3 base_v = base_plane.base2();
    Vector_3 top_u = top_plane.base1(); Vector_3 top_v = top_plane.base2();

    // Add center vertices
    typename Mesh::Vertex_index base_center_idx = cylinder.add_vertex(base_center);
    typename Mesh::Vertex_index top_center_idx = cylinder.add_vertex(top_center);

    std::vector<typename Mesh::Vertex_index> base_vertices;
    std::vector<typename Mesh::Vertex_index> top_vertices;
    base_vertices.reserve(num_segments);
    top_vertices.reserve(num_segments);

    // Generate vertices on the circles
    for (int i = 0; i < num_segments; ++i) {
        double angle = 2.0 * CGAL_PI * i / num_segments;
        FT cos_angle = FT(std::cos(angle)); FT sin_angle = FT(std::sin(angle));
        Point_3 base_pt = base_center + radius * (base_u * cos_angle + base_v * sin_angle);
        Point_3 top_pt = top_center + radius * (top_u * cos_angle + top_v * sin_angle);
        base_vertices.push_back(cylinder.add_vertex(base_pt));
        top_vertices.push_back(cylinder.add_vertex(top_pt));
    }

    // Create side faces - Explicitly use vertex 0 for the seam
    for (int i = 0; i < num_segments; ++i) {
        int next_i = (i + 1) % num_segments;
        cylinder.add_face(base_vertices[i], base_vertices[next_i], top_vertices[next_i]);
        cylinder.add_face(base_vertices[i], top_vertices[next_i], top_vertices[i]);
    }

    // Create base faces - Explicitly use vertex 0 for the seam
    for (int i = 0; i < num_segments; ++i) {
        int next_i = (i + 1) % num_segments;
        cylinder.add_face(base_center_idx, base_vertices[i], base_vertices[next_i]);
    }

    // Create top faces - Explicitly use vertex 0 for the seam
    for (int i = 0; i < num_segments; ++i) {
        int next_i = (i + 1) % num_segments;
        cylinder.add_face(top_center_idx, top_vertices[next_i], top_vertices[i]);
    }

    // Add Repair Steps INSIDE make_cylinder
    // 1. Remove degenerate faces first
    std::size_t removed_deg_faces = CGAL::Polygon_mesh_processing::remove_degenerate_faces(cylinder);
    if (removed_deg_faces > 0) {
        std::cout << "    make_cylinder: Removed " << removed_deg_faces << " degenerate faces." << std::endl;
    }

    // 2. Stitch borders (might help after removing degenerate faces)
    unsigned int borders_stitched = CGAL::Polygon_mesh_processing::stitch_borders(cylinder);
     if (borders_stitched > 0) {
        std::cout << "    make_cylinder: Stitched " << borders_stitched << " borders internally." << std::endl;
    }

    return cylinder;
}

// Add a helper function to load and normalize a pre-made cylinder STL
Mesh load_cylinder_stl(const std::string& stl_path) {
    Mesh cylinder;
    if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(stl_path, cylinder) || cylinder.is_empty()) {
        std::cerr << "Error: Failed to load cylinder STL from " << stl_path << std::endl;
        return Mesh();
    }

    std::cout << "Loaded cylinder STL successfully. Vertices: " << cylinder.number_of_vertices()
              << ", Faces: " << cylinder.number_of_faces() << std::endl;

    // Normalize the STL coordinates
    Iso_cuboid_3 bbox = CGAL::bounding_box(cylinder.points().begin(), cylinder.points().end());
    Point_3 min_bounds = bbox.min();
    Point_3 max_bounds = bbox.max();

    // Translate the STL so its base center is at the origin (0, 0, 0)
    // Assuming the cylinder is aligned with Z-axis in the STL
    FT center_x = (min_bounds.x() + max_bounds.x()) / FT(2);
    FT center_y = (min_bounds.y() + max_bounds.y()) / FT(2);
    Vector_3 translation_to_origin = Vector_3(-center_x, -center_y, -min_bounds.z());

    CGAL::Aff_transformation_3<Kernel> translate_to_origin(CGAL::TRANSLATION, translation_to_origin);
    CGAL::Polygon_mesh_processing::transform(translate_to_origin, cylinder);

    // Verify the new bounding box
    bbox = CGAL::bounding_box(cylinder.points().begin(), cylinder.points().end());
    std::cout << "Normalized STL bounding box: Min = " << bbox.min() << ", Max = " << bbox.max() << std::endl;

    // Ensure the cylinder is oriented along the Z-axis
    FT height = bbox.zmax() - bbox.zmin(); // Use bbox methods
    if (height <= FT(0)) {
        std::cerr << "Error: Loaded STL has invalid height after normalization." << std::endl;
        return Mesh();
    }

    return cylinder;
}


int main(int argc, char* argv[]) {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input,i", po::value<std::string>()->required(), "input STL file path")
        ("output,o", po::value<std::string>()->default_value("mold_output"), "output directory path")
        ("padding", po::value<double>()->default_value(6.0), "padding around the object (mm)")
        ("peg_radius", po::value<double>()->default_value(2.0), "radius of alignment pegs (mm)")
        ("peg_height", po::value<double>()->default_value(6.0), "height of alignment pegs (mm)")
        ("corner_offset", po::value<double>()->default_value(5.0), "offset of pegs/holes from corners (mm)")
    ;

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }
        po::notify(vm);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << desc << std::endl;
        return 1;
    }

    std::string input_file = vm["input"].as<std::string>();
    std::string output_dir = vm["output"].as<std::string>();
    FT padding = FT(vm["padding"].as<double>());
    FT peg_radius = FT(vm["peg_radius"].as<double>());
    FT peg_height = FT(vm["peg_height"].as<double>());
    FT corner_offset = FT(vm["corner_offset"].as<double>());

    // Create Output Directory if it doesn't exist
    try {
        fs::create_directories(output_dir);
        std::cout << "Using output directory: " << output_dir << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error creating output directory: " << e.what() << std::endl;
        return 1;
    }

    // Load Input Mesh for the part to be molded
    std::cout << "Input STL: " << input_file << std::endl;
    Mesh input_mesh;
    std::cout << "Loading mesh: " << input_file << "..." << std::endl;
    if (!CGALIO::read_polygon_mesh(input_file, input_mesh) || input_mesh.is_empty()) {
        std::cerr << "Error: Cannot read input file " << input_file << std::endl;
        return 1;
    }
    std::cout << "Mesh loaded successfully. Vertices: " << input_mesh.number_of_vertices()
              << ", Faces: " << input_mesh.number_of_faces() << std::endl;


    // Calculate Positions and Dimensions
    std::cout << "Calculating positions..." << std::endl;
    Iso_cuboid_3 original_bbox = CGAL::bounding_box(input_mesh.points().begin(), input_mesh.points().end());
    Point_3 original_min_bounds = original_bbox.min();
    Point_3 original_max_bounds = original_bbox.max();
    Point_3 original_center = CGAL::midpoint(original_min_bounds, original_max_bounds);
    Vector_3 mesh_size = original_max_bounds - original_min_bounds;

    std::cout << "  Original Bounds Min: " << original_min_bounds << std::endl;
    std::cout << "  Original Bounds Max: " << original_max_bounds << std::endl;
    std::cout << "  Original Center: " << original_center << std::endl;
    std::cout << "  Mesh Size: " << mesh_size << std::endl;

    Vector_3 box_dims = mesh_size + Vector_3(padding * FT(2), padding * FT(2), padding * FT(2));
    std::cout << "  Box dimensions with padding: " << box_dims << std::endl;

    // Translate Input Mesh
    Vector_3 translate_vector(
        -original_center.x(),
        -original_center.y(),
        padding - original_min_bounds.z() // Align bottom of original mesh with padding offset from box bottom
    );
    std::cout << "  Translation vector: " << translate_vector << std::endl;
    Aff_transformation_3 mesh_translation(CGAL::TRANSLATION, translate_vector);
    PMP::transform(mesh_translation, input_mesh);

    Iso_cuboid_3 translated_bbox = CGAL::bounding_box(input_mesh.points().begin(), input_mesh.points().end());
    std::cout << "  New bounds Min: " << translated_bbox.min() << std::endl;
    std::cout << "  New bounds Max: " << translated_bbox.max() << std::endl;
    std::cout << "  New Center: " << CGAL::midpoint(translated_bbox.min(), translated_bbox.max()) << std::endl;


    // Create Mold Box Halves, Split Plane, and Pegs
    FT split_z = padding + mesh_size.z() / FT(2); // Z-coordinate of the splitting plane
    std::cout << "Creating separate top and bottom mold halves with cavities..." << std::endl;
    std::cout << "  Split plane z position: " << split_z << std::endl;

    // Bottom box: from z=0 to z=split_z
    Point_3 bottom_box_center(0, 0, split_z / FT(2));
    Vector_3 bottom_box_dims(box_dims.x(), box_dims.y(), split_z);
    Mesh bottom_box = make_triangulated_cube(bottom_box_center, bottom_box_dims);
    std::cout << "  Bottom box center: " << bottom_box_center << ", dims: " << bottom_box_dims << std::endl;

    // Top box: from z=split_z to z=box_dims.z()
    FT top_box_height = box_dims.z() - split_z;
    Point_3 top_box_center(0, 0, split_z + top_box_height / FT(2));
    Vector_3 top_box_dims(box_dims.x(), box_dims.y(), top_box_height);
    Mesh top_box = make_triangulated_cube(top_box_center, top_box_dims);
    std::cout << "  Top box center: " << top_box_center << ", dims: " << top_box_dims << std::endl;

    // Ensure Orientations are Correct
    std::cout << "Checking and fixing mesh orientations..." << std::endl;
    if (!PMP::is_outward_oriented(top_box)) {
        std::cout << "  Fixing top box orientation..." << std::endl;
        PMP::orient(top_box);
    }
    if (!PMP::is_outward_oriented(bottom_box)) {
        std::cout << "  Fixing bottom box orientation..." << std::endl;
        PMP::orient(bottom_box);
    }
    if (!PMP::is_outward_oriented(input_mesh)) {
        std::cout << "  Fixing input mesh orientation..." << std::endl;
        PMP::orient(input_mesh);
    }

    // Create Bottom Half Mold 
    std::cout << "Creating bottom mold half..." << std::endl;
    Mesh bottom_half;
    try {
        std::cout << "  Performing boolean subtraction for bottom half..." << std::endl;
        Mesh box_copy_b = bottom_box;
        Mesh part_copy_b = input_mesh;
        bool bottom_ok = PMP::corefine_and_compute_difference(box_copy_b, part_copy_b, bottom_half);
        if (!bottom_ok || bottom_half.is_empty()) {
            std::cerr << "Error: Failed to create bottom mold half with cavity." << std::endl;
            return 1;
        } else {
            std::cout << "  Bottom half created successfully with " << bottom_half.number_of_vertices()
                      << " vertices and " << bottom_half.number_of_faces() << " faces" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Exception creating bottom half: " << e.what() << std::endl;
        return 1;
    }

    // Create Top Half Mold
    std::cout << "Creating top mold half..." << std::endl;
    Mesh top_half;
    try {
        std::cout << "  Performing boolean subtraction for top half..." << std::endl;
        Mesh box_copy_t = top_box;
        Mesh part_copy_t = input_mesh;
        bool top_ok = PMP::corefine_and_compute_difference(box_copy_t, part_copy_t, top_half);
        if (!top_ok || top_half.is_empty()) {
            std::cerr << "Error: Failed to create top mold half with cavity." << std::endl;
            return 1;
        } else {
            std::cout << "  Top half created successfully with " << top_half.number_of_vertices()
                      << " vertices and " << top_half.number_of_faces() << " faces" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Exception creating top half: " << e.what() << std::endl;
        return 1;
    }

    // Save Initial Mold Halves (Optional Debug)
    // CGALIO::write_STL(fs::path(output_dir) / "mold_top_initial_cpp.stl", top_half);
    // CGALIO::write_STL(fs::path(output_dir) / "mold_bottom_initial_cpp.stl", bottom_half);

    // Prepare Pegs and Holes, and Add to Mold Halves
    std::cout << "Adding alignment pegs and holes..." << std::endl;
    std::vector<Mesh> pegs;
    std::vector<Mesh> holes;
    std::string cylinder_stl_path = "../assets/peg.stl"; // Path relative to the build directory

    // Calculate positions for alignment pegs/holes
    FT box_width = box_dims.x();
    FT box_depth = box_dims.y();
    std::vector<Point_3> corner_positions = {
        Point_3(-box_width / 2 + corner_offset, -box_depth / 2 + corner_offset, split_z), // Front-left
        Point_3( box_width / 2 - corner_offset, -box_depth / 2 + corner_offset, split_z), // Front-right
        Point_3( box_width / 2 - corner_offset,  box_depth / 2 - corner_offset, split_z), // Back-right
        Point_3(-box_width / 2 + corner_offset,  box_depth / 2 - corner_offset, split_z)  // Back-left
    };

    // Load and position pegs
    int peg_idx_gen = 0;
    for (const auto& pos : corner_positions) {
        Mesh peg = load_cylinder_stl(cylinder_stl_path);
        if (peg.is_empty()) {
            std::cerr << "    ERROR: Failed to load cylinder STL for peg " << peg_idx_gen << std::endl;
            peg_idx_gen++;
            continue;
        }
        // Translate the normalized STL (base at origin) so its TOP aligns with pos (z=split_z)
        Vector_3 translation = pos - Point_3(0, 0, peg_height); // Translate origin so top is at pos
        CGAL::Aff_transformation_3<Kernel> translate(CGAL::TRANSLATION, translation);
        PMP::transform(translate, peg);

        if (!CGAL::is_closed(peg)) {
             std::cerr << "    ERROR: Peg " << peg_idx_gen << " STL is NOT closed after loading/transforming." << std::endl;
             peg_idx_gen++;
             continue;
        }
        if (!PMP::is_outward_oriented(peg)) PMP::orient(peg);

        // *** DEBUG: Save individual translated peg ***
        //std::string peg_filename = "debug_peg_" + std::to_string(peg_idx_gen) + ".stl";
        //CGAL::IO::write_STL(fs::path(output_dir) / peg_filename, peg);
        

        pegs.push_back(peg);
        peg_idx_gen++;
    }

    // Load and position holes
    int hole_idx_gen = 0;
    for (const auto& pos : corner_positions) {
        Mesh hole = load_cylinder_stl(cylinder_stl_path);
        if (hole.is_empty()) {
            std::cerr << "    ERROR: Failed to load cylinder STL for hole " << hole_idx_gen << std::endl;
            hole_idx_gen++;
            continue;
        }
        // Translate the normalized STL (base at origin) so its *top* aligns with the hole position (split_z)
        // The STL height is peg_height
        Vector_3 translation = pos - Point_3(0, 0, peg_height); // Translate origin so top is at pos
        CGAL::Aff_transformation_3<Kernel> translate(CGAL::TRANSLATION, translation);
        PMP::transform(translate, hole);

        if (!CGAL::is_closed(hole)) {
             std::cerr << "    ERROR: Hole " << hole_idx_gen << " STL is NOT closed after loading/transforming." << std::endl;
             hole_idx_gen++;
             continue;
        }
        if (!PMP::is_outward_oriented(hole)) PMP::orient(hole);
        holes.push_back(hole);
        hole_idx_gen++;
    }

    std::cout << "  Generated " << pegs.size() << " pegs and " << holes.size() << " holes." << std::endl;


    // Perform Boolean Operations 
    Mesh top_with_pegs = top_half; // Start with the ORIGINAL top half
    Mesh bottom_with_holes = bottom_half; // Start with the original bottom half

    // Add pegs to top half
    std::cout << "  Adding alignment pegs to top half..." << std::endl;
    int peg_idx_op = 0;
    bool any_peg_failed = false; // Flag to track failures
    for (auto& peg : pegs) { // Pegs are now translated correctly for the UNROTATED top_half
        Mesh temp_result_p;
        std::cout << "    Attempting union with peg " << peg_idx_op << "..." << std::endl;
        try {
            bool union_ok = PMP::corefine_and_compute_union(top_with_pegs, peg, temp_result_p);
            if (union_ok && !temp_result_p.is_empty() && CGAL::is_valid_polygon_mesh(temp_result_p)) {
                top_with_pegs = temp_result_p; // Update the intermediate result
                std::cout << "      Successfully added peg " << peg_idx_op << std::endl;
            } else {
                // Enhanced warning
                std::cerr << "      Warning: Failed to add peg " << peg_idx_op
                          << " (union_ok=" << union_ok
                          << ", empty=" << temp_result_p.is_empty()
                          << ", valid=" << (!temp_result_p.is_empty() && CGAL::is_valid_polygon_mesh(temp_result_p))
                          << ")" << std::endl;
                any_peg_failed = true;
            }
        } catch (const std::exception& e) {
            std::cerr << "      Exception adding peg " << peg_idx_op << ": " << e.what() << std::endl;
            any_peg_failed = true;
        }
        peg_idx_op++;
    }

    // Add holes to bottom half
    std::cout << "  Adding alignment holes to bottom half..." << std::endl;
    int hole_idx_op = 0;
    for (auto& hole : holes) {
        Mesh temp_result_h;
        std::cout << "    Attempting difference with hole " << hole_idx_op << "..." << std::endl;
        try {
            bool diff_ok = PMP::corefine_and_compute_difference(bottom_with_holes, hole, temp_result_h);
            if (diff_ok && !temp_result_h.is_empty() && CGAL::is_valid_polygon_mesh(temp_result_h)) {
                bottom_with_holes = temp_result_h;
                std::cout << "      Successfully added hole " << hole_idx_op << std::endl;
            } else {
                std::cerr << "      Warning: Failed to add hole " << hole_idx_op << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "      Exception adding hole " << hole_idx_op << ": " << e.what() << std::endl;
        }
        hole_idx_op++;
    }

    std::cout << "Alignment features processing finished." << std::endl;

    // *** DEBUG: Save top_with_pegs BEFORE rotation ***
    //if (!any_peg_failed) { // Only save if unions seemed okay
    //     CGAL::IO::write_STL(fs::path(output_dir) / "debug_top_with_pegs_pre_rotation.stl", top_with_pegs);
    //}


    // Rotate the Top Mold Half After Adding Pegs
    std::cout << "Rotating top mold half with pegs by 180 degrees..." << std::endl;
    CGAL::Aff_transformation_3<Kernel> rotate_180_around_X(
        FT(1), FT(0),  FT(0), FT(0),
        FT(0), FT(-1), FT(0), FT(0),
        FT(0), FT(0),  FT(-1), FT(0)
    );
    PMP::transform(rotate_180_around_X, top_with_pegs); // Rotate the combined mesh

    // Ensure orientation after rotation
    if (!PMP::is_outward_oriented(top_with_pegs)) {
        std::cout << "  Fixing top mold half with pegs orientation after rotation..." << std::endl;
        PMP::orient(top_with_pegs);
    }

    // TODO: Auto-pacement. Here We Will Add The Sprue and Vents
    // This requires some additional logic to determine the best locations
    // and sizes for the sprue and vents based on the geometry of the part.
    // This involves finding a high point on the part and creating a channel
    // for the sprue, applying a boolean operation to add it to the top half
    // using a nozzle or similar shape defined in assets/nozzle.stl.
    // Vents can be cylindrical holes or channels that allow air to escape,
    // positioned where the part is lowest or where air might be trapped.

    // Create Sprue and Vents
    // This requires some additional logic to determine the best locations
    // and sizes for the sprue and vents based on the geometry of the part.
    // This involves finding a high point on the part and creating a channel
    // for the sprue, applying a boolean operation to add it to the top half
    // using a nozzle or similar shape defined in assets/nozzle.stl.
    // Vents can be cylindrical holes or channels that allow air to escape,
    // positioned where the part is lowest or where air might be trapped.
    //
    // Something like this for sprue:
    // FT min_z = +∞; Point_3 pv;
    // for (auto p : input_mesh.points())
    //   if (p.z() < min_z) { min_z = p.z(); pv = p; }
    // sprue_center_xy = (pv.x(), pv.y()), vent_top_z = split_z + lid_thickness
    //
    // Then for vents:
    // find roof plane z-height 
    // FT roofZ = −∞;
    // for (auto& p : input_mesh.points())
    //   roofZ = max(roofZ, p.z());
    //
    // Find depth:
    // vector<pair<FT,Point_3>> depths;
    // for (auto& p : mesh.points())
    //   depths.emplace_back(roofZ - p.z(), p);
    //
    // Sort by depth descending:
    // sort(depths.begin(), depths.end(),
    //    [](auto&a,auto&b){ return a.first>b.first; });
    //
    // Then depths[0..K] are the deepest pockets, for small parts
    // probably only 1 vent is needed, so this can be sized by part
    // but may take some physical experimentation to get the ratio
    // and also depends on the material.
    // END TODO


    // Save Final Mold Halves
    // Assign the results of boolean operations (and rotation) back before saving
    top_half = top_with_pegs; // Assign the rotated mesh with pegs
    bottom_half = bottom_with_holes;

    std::cout << "Saving final mold halves..." << std::endl;
    if (!CGAL::IO::write_STL(fs::path(output_dir) / "mold_top.stl", top_half)) {
        std::cerr << "Error writing final top mold half." << std::endl;
    } else {
        std::cout << "  Final top half saved successfully." << std::endl;
    }

    if (!CGAL::IO::write_STL(fs::path(output_dir) / "mold_bottom.stl", bottom_half)) {
        std::cerr << "Error writing final bottom mold half." << std::endl;
    } else {
        std::cout << "  Final bottom half saved successfully." << std::endl;
    }

    std::cout << "Processing finished." << std::endl;
    return 0;
}
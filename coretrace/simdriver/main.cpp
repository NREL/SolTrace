/**
 * @file main.cpp
 * @brief Command-line driver for SolTrace ray tracing.
 *
 * Reads a JSON or .stinput file to configure SimulationData, runs the ray tracer,
 * and writes the ray interaction records to a CSV file.
 *
 * Usage:
 *   simdriver <input.json|input.stinput> [<output.csv>] [options]
 *
 * Options:
 *   --threads <n>   Number of parallel threads (default: 1)
 *   --rays <n>      Override the number of rays from the JSON file
 *   --no-output     Skip result retrieval and CSV output (output file not
 *                   required when this flag is set)
 *   --no-csv        Retrieve results but skip writing the CSV file (output
 *                   file argument not required when this flag is set)
 *   --embree        Use the Embree runner (only available if built with
 *                   SOLTRACE_BUILD_EMBREE_SUPPORT=ON; falls back to native
 *                   runner with a warning if Embree support is absent)
 *   --optix         Use the OptiX runner (only available if built with
 *                   SOLTRACE_BUILD_OPTIX_SUPPORT=ON; falls back to native
 *                   runner with a warning if OptiX support is absent)
 *   --verbose       Enable verbose logging in the OptiX runner
 */

#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <string>

#include "native_runner.hpp"
#include "simulation_data_export.hpp"
#include "simulation_result_export.hpp"

#ifdef SOLTRACE_EMBREE_SUPPORT
#include "embree_runner.hpp"
#endif

#ifdef SOLTRACE_OPTIX_SUPPORT
#include "optix_runner.hpp"
#endif

// using SolTrace::Data::SimulationData;
// using SolTrace::Data::SimulationParameters;
using SolTrace::NativeRunner::NativeRunner;
using SolTrace::Result::SimulationResult;
using SolTrace::Runner::RunnerStatus;

static void print_usage(const char *prog)
{
    std::cerr
        << "Usage: " << prog << " <input.json|input.stinput> [<output_filename>] [options]\n"
        << "Please include the .json extension for the input.json argument, while"
        << " excluding a file extension for the output filename. The output file "
        << "extension will be determined based on the level specified (see below).\n"
        << "\n"
        << "Options:\n"
        << "  --threads <n>   Number of threads (default: 1)\n"
        << "  --rays <n>      Override number of rays specified in the JSON file\n"
        << "  --no-output     Skip result retrieval and CSV output\n"
        << "                  (output file argument not required with this flag)\n"
        << "  --no-csv        Retrieve results but skip writing the CSV file\n"
        << "                  (output file argument not required with this flag)\n"
        << "  --level <n>     Runner reporting level (default: 0, see SolTrace::Runner::RunnerStatistics for available levels)\n"
#ifdef SOLTRACE_EMBREE_SUPPORT
        << "  --embree        Use Embree runner instead of the native runner\n"
        << "                  (requires SOLTRACE_BUILD_EMBREE_SUPPORT=ON at build time)\n"
#endif
#ifdef SOLTRACE_OPTIX_SUPPORT
        << "  --optix         Use OptiX runner instead of the native runner\n"
        << "                  (requires SOLTRACE_BUILD_OPTIX_SUPPORT=ON at build time)\n"
#endif
        << "  -v, --verbose   Enable verbose logging in the OptiX runner\n"
        << "  -t, --timing    Enable timing logging in the OptiX runner\n"
        ;
}

int main(int argc, char *argv[])
{
    using SolTrace::Runner::RunnerStatistics;
    // double check that this works with no-output/no-csv flags
    if (argc < 3)
    {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    // Pre-scan for --no-output and --no-csv so we know whether output_file is required
    bool skip_output = false;
    bool skip_csv = false;
    for (int i = 2; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--no-output") skip_output = true;
        else if (a == "--no-csv") skip_csv = true;
    }

    const bool file_optional = skip_output || skip_csv;

    const std::string input_file = argv[1];

    // argv[2], if present and not a flag (does not start with --), is treated as
    // the output file path.  This allows the user to supply an output path even
    // when --no-output or --no-csv is also present without it being mis-parsed
    // as an unknown option.
    const bool has_output_arg = (argc >= 3) && (std::string(argv[2]).rfind("--", 0) != 0);
    const std::string output_file = has_output_arg ? argv[2] : "";
    const int opts_start = has_output_arg ? 3 : 2;

    if (!file_optional && !has_output_arg)
    {
        std::cerr << "Error: output file is required unless --no-output or --no-csv is specified\n";
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    int num_threads = 1;
    long long num_rays_override = -1; // -1 means use what the JSON specifies
    int level = 0;
    bool use_embree = false;
    bool use_optix = false;
    bool verbose = false;
    bool log_timing = false;

    for (int i = opts_start; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--threads")
        {
            if (i + 1 >= argc)
            {
                std::cerr << "Error: --threads requires an argument\n";
                return EXIT_FAILURE;
            }
            try
            {
                num_threads = std::stoi(argv[++i]);
            }
            catch (...)
            {
                std::cerr << "Error: invalid thread count '" << argv[i] << "'\n";
                return EXIT_FAILURE;
            }
            if (num_threads < 1)
            {
                std::cerr << "Error: thread count must be >= 1\n";
                return EXIT_FAILURE;
            }
        }
        else if (arg == "--rays")
        {
            if (i + 1 >= argc)
            {
                std::cerr << "Error: --rays requires an argument\n";
                return EXIT_FAILURE;
            }
            try
            {
                num_rays_override = std::stoll(argv[++i]);
            }
            catch (...)
            {
                std::cerr << "Error: invalid ray count '" << argv[i] << "'\n";
                return EXIT_FAILURE;
            }
            if (num_rays_override < 1)
            {
                std::cerr << "Error: ray count must be >= 1\n";
                return EXIT_FAILURE;
            }
        }
        else if (arg == "--no-output" || arg == "--no-csv")
        {
            // already handled in pre-scan; skip here
        }
        else if (arg == "--level")
        {
            if (i + 1 >= argc)
            {
                std::cerr << "Error: --level requires an argument\n";
                return EXIT_FAILURE;
            }
            try
            {
                level = std::stoi(argv[++i]);
            }
            catch (...)
            {
                std::cerr << "Error: invalid level '" << argv[i] << "'\n";
                return EXIT_FAILURE;
            }
            if (level < RunnerStatistics::RAY_RECORDS || level > RunnerStatistics::ALL)
            {
                std::cerr << "Error: level must map to SolTrace::Runner::RunnerStatistics\n";
                return EXIT_FAILURE;
            }
        }
        else if (arg == "--embree")
        {
            use_embree = true;
        }
        else if (arg == "--optix")
        {
            use_optix = true;
        }
        else if (arg == "-v" || arg == "--verbose")
        {
            verbose = true;
        }
        else if (arg == "-t" || arg == "--timing")
        {
            log_timing = true;
        }
        else
        {
            std::cerr << "Error: unknown option '" << arg << "'\n";
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    // -------------------------------------------------------------------------
    // Load simulation data from JSON or .stinput file
    // -------------------------------------------------------------------------
    SimulationData simData;
    {
        // Determine format by extension
        auto ends_with = [](const std::string &s, const std::string &suffix)
        {
            return s.size() >= suffix.size() &&
                   s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
        };
        const bool is_stinput = ends_with(input_file, ".stinput");
        const bool is_json = ends_with(input_file, ".json");

        if (!is_stinput && !is_json)
        {
            std::cerr << "Error: unrecognised input file extension (expected .json or .stinput): "
                      << input_file << "\n";
            return EXIT_FAILURE;
        }

        std::cout << "Loading simulation data from: " << input_file << "...\n";
        auto t_load_start = std::chrono::steady_clock::now();

        if (is_json)
        {
            try
            {
                simData.import_json_file(input_file);
            }
            catch (const std::exception &e)
            {
                std::cerr << "Error loading JSON file: " << e.what() << "\n";
                return EXIT_FAILURE;
            }
        }
        else // .stinput
        {
            if (!simData.import_from_file(input_file))
            {
                std::cerr << "Error loading .stinput file: " << input_file << "\n";
                return EXIT_FAILURE;
            }
        }

        auto t_load_end = std::chrono::steady_clock::now();
        std::cout << "  Loaded in "
                  << std::chrono::duration<double>(t_load_end - t_load_start).count()
                  << " s\n"
                  << "  Elements loaded: " << simData.get_number_of_elements() << "\n";
    }

    // Override ray count if the user requested it
    if (num_rays_override > 0)
    {
        SimulationParameters &params = simData.get_simulation_parameters();
        params.number_of_rays = static_cast<uint_fast64_t>(num_rays_override);
        params.max_number_of_rays = params.number_of_rays * 100;
        std::cout << "Overriding ray count to " << params.number_of_rays << "\n";
    }

    // -------------------------------------------------------------------------
    // Set up and run the simulation
    // -------------------------------------------------------------------------
    RunnerStatus sts;
    SimulationResult result;

#ifdef SOLTRACE_EMBREE_SUPPORT
    if (use_embree)
    {
        SolTrace::EmbreeRunner::EmbreeRunner runner;

        sts = runner.initialize();
        if (sts != RunnerStatus::SUCCESS)
        {
            std::cerr << "Error: failed to initialize Embree runner\n";
            return EXIT_FAILURE;
        }

        runner.set_number_of_threads(static_cast<uint_fast64_t>(num_threads));
        std::cout << "Using Embree runner with " << num_threads << " thread(s)\n";

        std::cout << "Setting up simulation...\n";
        auto t_setup_start = std::chrono::steady_clock::now();
        runner.disable_stages();
        sts = runner.setup_simulation(&simData);
        auto t_setup_end = std::chrono::steady_clock::now();
        if (sts != RunnerStatus::SUCCESS)
        {
            std::cerr << "Error: Embree runner setup failed\n";
            return EXIT_FAILURE;
        }
        std::cout << "  Setup in "
                  << std::chrono::duration<double>(t_setup_end - t_setup_start).count()
                  << " s\n";

        std::cout << "Running simulation...\n";
        auto t_run_start = std::chrono::steady_clock::now();
        sts = runner.run_simulation();
        auto t_run_end = std::chrono::steady_clock::now();
        if (sts != RunnerStatus::SUCCESS)
        {
            std::cerr << "Error: simulation failed\n";
            return EXIT_FAILURE;
        }
        std::cout << "  Completed in "
                  << std::chrono::duration<double>(t_run_end - t_run_start).count()
                  << " s\n";
        std::cout << "  Rays launched: " << runner.get_number_rays_launched() << "\n";
        std::cout << "  Rays traced:   " << runner.get_number_rays_traced() << "\n";

        if (!skip_output)
        {
            std::cout << "Retrieving results...\n";
            auto t_report_start = std::chrono::steady_clock::now();
            sts = runner.report_simulation(&result, level);
            auto t_report_end = std::chrono::steady_clock::now();
            if (sts != RunnerStatus::SUCCESS)
            {
                std::cerr << "Error: failed to collect simulation results\n";
                return EXIT_FAILURE;
            }
            std::cout << "  Retrieved in "
                      << std::chrono::duration<double>(t_report_end - t_report_start).count()
                      << " s\n";
        }
        else
        {
            std::cout << "Skipping result retrieval (--no-output).\n";
        }
    }
    else
#endif
#ifdef SOLTRACE_OPTIX_SUPPORT
        if (use_optix)
    {
        OptixRunner runner;

        sts = runner.initialize();
        if (sts != RunnerStatus::SUCCESS)
        {
            std::cerr << "Error: failed to initialize OptiX runner\n";
            return EXIT_FAILURE;
        }

        runner.set_verbose(verbose);

        std::cout << "Using OptiX runner\n";

        std::cout << "Setting up simulation...\n";
        auto t_setup_start = std::chrono::steady_clock::now();
        sts = runner.setup_simulation(&simData);
        auto t_setup_end = std::chrono::steady_clock::now();
        if (sts != RunnerStatus::SUCCESS)
        {
            std::cerr << "Error: OptiX runner setup failed\n";
            return EXIT_FAILURE;
        }
        std::cout << "  Setup in "
                  << std::chrono::duration<double>(t_setup_end - t_setup_start).count()
                  << " s\n";

        std::cout << "Running simulation...\n";
        auto t_run_start = std::chrono::steady_clock::now();
        sts = runner.run_simulation();
        auto t_run_end = std::chrono::steady_clock::now();
        if (sts != RunnerStatus::SUCCESS)
        {
            std::cerr << "Error: simulation failed\n";
            return EXIT_FAILURE;
        }
        std::cout << "  Completed in "
                  << std::chrono::duration<double>(t_run_end - t_run_start).count()
                  << " s\n";
        std::cout << "  Rays launched: " << runner.get_number_rays_launched() << "\n";
        std::cout << "  Rays traced:   " << runner.get_number_rays_traced() << "\n";

        if (!skip_output)
        {
            std::cout << "Retrieving results at level " << level << "...\n";
            auto t_report_start = std::chrono::steady_clock::now();
            sts = runner.report_simulation(&result, level);
            auto t_report_end = std::chrono::steady_clock::now();
            if (sts != RunnerStatus::SUCCESS)
            {
                std::cerr << "Error: failed to collect simulation results\n";
                return EXIT_FAILURE;
            }
            std::cout << "  Retrieved in "
                      << std::chrono::duration<double>(t_report_end - t_report_start).count()
                      << " s\n";
        }
        else
        {
            std::cout << "Skipping result retrieval (--no-output).\n";
        }

        if (verbose || log_timing)
        {
            runner.print_timing();
        }
    }
    else
#endif
    {
#ifndef SOLTRACE_EMBREE_SUPPORT
        if (use_embree)
        {
            std::cerr << "Warning: this build does not include Embree support. "
                         "Falling back to the native runner.\n";
        }
#endif
#ifndef SOLTRACE_OPTIX_SUPPORT
        if (use_optix)
        {
            std::cerr << "Warning: this build does not include OptiX support. "
                         "Falling back to the native runner.\n";
        }
#endif

        NativeRunner runner;

        sts = runner.initialize();
        if (sts != RunnerStatus::SUCCESS)
        {
            std::cerr << "Error: failed to initialize native runner\n";
            return EXIT_FAILURE;
        }

        runner.set_number_of_threads(static_cast<uint_fast64_t>(num_threads));
        std::cout << "Using native runner with " << num_threads << " thread(s)\n";

        std::cout << "Setting up simulation...\n";
        auto t_setup_start = std::chrono::steady_clock::now();
        sts = runner.setup_simulation(&simData);
        auto t_setup_end = std::chrono::steady_clock::now();
        if (sts != RunnerStatus::SUCCESS)
        {
            std::cerr << "Error: native runner setup failed\n";
            return EXIT_FAILURE;
        }
        std::cout << "  Setup in "
                  << std::chrono::duration<double>(t_setup_end - t_setup_start).count()
                  << " s\n";

        std::cout << "Running simulation...\n";
        auto t_run_start = std::chrono::steady_clock::now();
        sts = runner.run_simulation();
        auto t_run_end = std::chrono::steady_clock::now();
        if (sts != RunnerStatus::SUCCESS)
        {
            std::cerr << "Error: simulation failed\n";
            return EXIT_FAILURE;
        }
        std::cout << "  Completed in "
                  << std::chrono::duration<double>(t_run_end - t_run_start).count()
                  << " s\n";
        std::cout << "  Rays launched: " << runner.get_number_rays_launched() << "\n";
        std::cout << "  Rays traced:   " << runner.get_number_rays_traced() << "\n";

        if (!skip_output)
        {
            std::cout << "Retrieving results...\n";
            auto t_report_start = std::chrono::steady_clock::now();
            sts = runner.report_simulation(&result, level);
            auto t_report_end = std::chrono::steady_clock::now();
            if (sts != RunnerStatus::SUCCESS)
            {
                std::cerr << "Error: failed to collect simulation results\n";
                return EXIT_FAILURE;
            }
            std::cout << "  Retrieved in "
                      << std::chrono::duration<double>(t_report_end - t_report_start).count()
                      << " s\n";
        }
        else
        {
            std::cout << "Skipping result retrieval (--no-output).\n";
        }
    }

    // -------------------------------------------------------------------------
    // Write results to CSV / JSON
    // -------------------------------------------------------------------------
    if (!skip_csv && (level == RunnerStatistics::RAY_RECORDS || level == RunnerStatistics::ALL))
    {
        std::cout << "Writing " << result.get_number_of_records()
                  << " ray records to: " << output_file << ".csv ...\n";
        try
        {
            auto t_write_start = std::chrono::steady_clock::now();
            result.write_csv_file(output_file + ".csv");
            auto t_write_end = std::chrono::steady_clock::now();
            std::cout << "  Written in "
                      << std::chrono::duration<double>(t_write_end - t_write_start).count()
                      << " s\n";
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error writing CSV file: " << e.what() << "\n";
            return EXIT_FAILURE;
        }
    }
    else if (skip_csv)
    {
        std::cout << "Skipping CSV output (--no-csv).\n";
    }

    if (!skip_output && (level == RunnerStatistics::GROUPED_COUNTS || level == RunnerStatistics::ALL))
    {
        std::cout << "Writing " << result.get_number_of_groups()
            << " group results to: " << output_file << ".json ...\n";
        try
        {
            auto t_write_start = std::chrono::steady_clock::now();
            result.write_group_json_file(output_file + ".json");
            auto t_write_end = std::chrono::steady_clock::now();
            std::cout << "  Written in "
            << std::chrono::duration<double>(t_write_end - t_write_start).count()
            << " s\n";
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error writing JSON file: " << e.what() << "\n";
            return EXIT_FAILURE;
        }
    }
    else if (skip_output)
    {
        std::cout << "Skipping CSV output (--no-output).\n";
    }
    
    std::cout << "Done.\n";
    return EXIT_SUCCESS;
}

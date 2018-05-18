

#include "program_arguments.hpp"

int main(int argc, char* argv[])
{

    ////////////////////////////////////////////////////////////////////////////
    // PROCESS PROGRAM ARGUMENTS
    ////////////////////////////////////////////////////////////////////////////

    ProgramArguments pa;
    pa.Add("help", "--help", "false");
    pa.Add("input_file", "--input-file", "of_data.txt");
    pa.Print();
    pa.Process(argc, argv);

    std::string arg_input_file{pa.Get("input_file")};

    std::cout << "input_file=" << arg_input_file << std::endl;

    return 0;
}

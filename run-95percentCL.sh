#!/bin/bash
./build/main --fit-subrange false --energy-cut false --epsilon 0.06 > ./95_Percent_CL/FULLRANGE_NOCUT/log.txt
mv ./el_energy_log_dir/* ./95_Percent_CL/FULLRANGE_NOCUT/
mv ./el_energy_sum_log_dir/* ./95_Percent_CL/FULLRANGE_NOCUT/

./build/main --fit-subrange false --energy-cut true --epsilon 0.10 > ./95_Percent_CL/FULLRANGE_CUT/log.txt
mv ./el_energy_log_dir/* ./95_Percent_CL/FULLRANGE_CUT/
mv ./el_energy_sum_log_dir/* ./95_Percent_CL/FULLRANGE_CUT/

./build/main --fit-subrange true --energy-cut false --epsilon 0.33 > ./95_Percent_CL/SUBRANGE_NOCUT/log.txt
mv ./el_energy_log_dir/* ./95_Percent_CL/SUBRANGE_NOCUT/
mv ./el_energy_sum_log_dir/* ./95_Percent_CL/SUBRANGE_NOCUT/

./build/main --fit-subrange true --energy-cut true --epsilon 0.67 > ./95_Percent_CL/SUBRANGE_CUT/log.txt
mv ./el_energy_log_dir/* ./95_Percent_CL/SUBRANGE_CUT/
mv ./el_energy_sum_log_dir/* ./95_Percent_CL/SUBRANGE_CUT/

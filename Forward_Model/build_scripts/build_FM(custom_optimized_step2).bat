g++ -O3 -flto -fbranch-probabilities -Wall ../source/strand.cpp ../source/loop.cpp ../source/ion.cpp ../source/instrument.cpp ../source/forward.cpp ../../Radiation_Model/source/element.cpp ../../Radiation_Model/source/radiation.cpp ../../Resources/source/file.cpp ../../Resources/source/fitpoly.cpp ../source/main.cpp -o ../../Forward_Model.exe
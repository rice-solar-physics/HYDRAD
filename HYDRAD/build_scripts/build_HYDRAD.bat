g++ -O3 -flto -Wall ../source/main.cpp ../source/cell.cpp ../source/mesh.cpp ../source/eqns.cpp ../../Kinetic_Model/source/kinetic.cpp ../../Kinetic_Model/source/gamma.cpp ../../Heating_Model/source/heat.cpp ../../Radiation_Model/source/ionfrac.cpp ../../Radiation_Model/source/element.cpp ../../Radiation_Model/source/radiation.cpp ../../Radiation_Model/source/OpticallyThick/OpticallyThickIon.cpp ../../Radiation_Model/source/OpticallyThick/RadiativeRates.cpp ../../Resources/source/gammabeta.cpp ../../Resources/source/fitpoly.cpp ../../Resources/source/file.cpp -o ../../HYDRAD.exe
#ifndef _PARALLEL_HPP_
#define _PARALLEL_HPP_

#include "comobject.hpp"

int Separate_Read(std::string name, std::istringstream& is);
int Separate_Write(std::string name, std::ostringstream& os);

int Shared_Read(std::string name, std::istringstream& is);
int Shared_Write(std::string name, std::ostringstream& os);

#endif

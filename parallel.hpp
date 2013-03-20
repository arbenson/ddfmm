#ifndef _PARALLEL_HPP_
#define _PARALLEL_HPP_

#include "comobject.hpp"

using std::istringstream;
using std::ostringstream;

int Separate_Read(string name, istringstream& is);
int Separate_Write(string name, ostringstream& os);

int Shared_Read(string name, istringstream& is);
int Shared_Write(string name, ostringstream& os);

#endif

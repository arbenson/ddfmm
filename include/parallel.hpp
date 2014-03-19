/* Distributed Directional Fast Multipole Method
   Copyright (C) 2014 Austin Benson, Lexing Ying, and Jack Poulson

 This file is part of DDFMM.

    DDFMM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DDFMM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DDFMM.  If not, see <http://www.gnu.org/licenses/>. */
#ifndef _PARALLEL_HPP_
#define _PARALLEL_HPP_

#include "comobject.hpp"

int SeparateRead(std::string name, std::istringstream& is);
int SeparateWrite(std::string name, std::ostringstream& os);

int SharedRead(std::string name, std::istringstream& is);
int SharedWrite(std::string name, std::ostringstream& os);

#endif

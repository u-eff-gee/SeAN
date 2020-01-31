/*    
    This file is part of SeAN.

    SeAN is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SeAN is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SeAN.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once 

#include <chrono>
#include <iostream>
#include <string>

using std::chrono::high_resolution_clock;
using std::cout;
using std::endl;
using std::string;

class Status{

public:
    Status(){};
    ~Status(){};

    void start_timer(){ t_start = high_resolution_clock::now(); };

    static string status_prefix(){
        return "> STATUS ";
    }

protected:
    static high_resolution_clock::time_point t_start;

};
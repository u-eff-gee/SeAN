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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using std::chrono::high_resolution_clock;
using std::cout;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::endl;
using std::scientific;
using std::setprecision;
using std::setw;
using std::string;
using std::stringstream;

class Status{

public:
    Status() = default;
    ~Status() = default;

    static void start_timer(){ t_start = high_resolution_clock::now(); };

    static void print(const string status_message, const bool print_elapsed_time) {
        cout << status_prefix(print_elapsed_time) << status_message << endl;
    }

protected:
    static high_resolution_clock::time_point t_start;

    static string status_prefix(const bool print_elapsed_time){
        duration<double> delta_t = duration_cast<duration<double>>(high_resolution_clock::now() - t_start);
        stringstream sta_pre;
        sta_pre << "> STATUS ";
        if(print_elapsed_time){
            sta_pre << "( " << scientific << setprecision(4) << delta_t.count() << " s )";
        }
        sta_pre << " : ";
        return sta_pre.str();
    }

};
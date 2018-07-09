/*
 * threadguard.cpp
 *
 *  Created on: Mar 23, 2018
 *      Author: ex56yseb
 */

#include "thread_guard.h"

thread_guard::thread_guard(std::thread& t_): t(std::move(t_)){
}


thread_guard::~thread_guard(){
	if (t.joinable()) t.join();
}

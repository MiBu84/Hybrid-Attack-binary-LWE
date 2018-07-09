/*
 * threadguard.h
 *
 *  Created on: Mar 23, 2018
 *      Author: ex56yseb
 */

#ifndef INCLUDE_THREAD_GUARD_H_
#define INCLUDE_THREAD_GUARD_H_

#include <thread>

class thread_guard {
private:
	std::thread t;
public:
	thread_guard(std::thread&);
	~thread_guard();
};

#endif /* INCLUDE_THREAD_GUARD_H_ */

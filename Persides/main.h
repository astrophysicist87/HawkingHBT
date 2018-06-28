#ifndef MAIN_H
#define MAIN_H

extern const int smax;

inline int indexer_xinm(int n, int m)
{
	//range of t is 0<=t<=s ==> just use smax
	return ( n * smax + m );
}

#endif

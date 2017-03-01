#ifndef _UTILS_H_
#define _UTILS_H_

int idx(int3 c, int3 s)
{
	return (c.z * s.y + c.y) * s.x + c.x;
}

int4 idx4(int3 c, int3 s)
{
	int4 i;
	i.x = (c.z * s.y + c.y) * s.x + c.x; // (c.x, c.y, c.z)
	i.y = i.x + s.x;                     // (c.x, c.y + 1, c.z)
	i.z = i.x + s.x * s.y;               // (c.x, c.y, c.z + 1)
	i.w = i.z + s.x;                     // (c.x, c.y + 1, c.z + 1)
	return i;
}
#endif
#ifndef _UTILS_H_
#define _UTILS_H_

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#define PRECISION 0.0001

struct WorkItemInfo
{
	int3 local_size;
	int3 global_size;
	int3 group_id;
	int3 cell_id;
	int3 global_cell_id;
};

void initWorkItemInfo(struct WorkItemInfo *wi)
{
	wi->local_size  = (int3)(get_local_size(0), get_local_size(1), get_local_size(2));
	wi->global_size = (int3)(get_global_size(0), get_global_size(1), get_global_size(2));
	wi->group_id    = (int3)(get_group_id(0), get_group_id(1), get_group_id(2));
	wi->cell_id     = (int3)(get_local_id(0), get_local_id(1), get_local_id(2));
	wi->global_cell_id = wi->group_id * wi->local_size + wi->cell_id;
}

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
#include "nauty.h"
#include "naugroup.h"

DYNALLSTAT(int, allp, allp_sz);
DYNALLSTAT(int, id, id_sz);

int k;

void store_perm(int *p, int n, int *target_array)
{
  int i;
  for (i = 0; i < n; ++i) target_array[k*n + i] = p[i];
  ++k;
}

static void store_group_elts(
  levelrec *lr,
  int n,
  int level,
  int *target_array,
  int *before,
  int *after,
  int *id)
{
  int i, j, orbsize;
  int *p, *cr;
  cosetrec *coset;

  coset = lr[level].replist;
  orbsize = lr[level].orbitsize;

  for (j = 0; j < orbsize; ++j)
  {
    cr = (coset[j].rep == NULL ? NULL : coset[j].rep->p);
    if (before == NULL)
      p = cr;
    else if (cr == NULL)
      p = before;
    else
    {
      p = after;
      for (i = 0; i < n; ++i) p[i] = cr[before[i]];
    }

    if (level == 0)
    {
      store_perm((p == NULL ? id : p), n, target_array);
    }
    else
    {
      store_group_elts(lr, n, level-1, target_array, p, after+n, id);
    }
  }
}

void store_group(grouprec *grp,int *target_array)
{
  int i, depth, n;

  depth = grp->depth;
  n = grp->n;


  DYNALLOC1(int, id, id_sz, n, "malloc");
  for (i = 0; i < n; ++i) id[i] = i;

  if (depth == 0)
  {
    store_perm(id, n, target_array);
    return;
  }

  DYNALLOC1(int, allp, allp_sz, n*depth, "malloc");

  store_group_elts(grp->levelinfo, n, depth-1, target_array, NULL, allp, id);
  
}

void find_automorphism_group(graph *g, int n, int *target_array, int target_array_size)
{  
  k = 0;
  
  DYNALLSTAT(int, lab, lab_sz);
  DYNALLSTAT(int, ptn, ptn_sz);
  DYNALLSTAT(int, orbits, orbits_sz);
  static DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;

  int m, v;
  grouprec *group;

  options.userautomproc = groupautomproc;
  options.userlevelproc = grouplevelproc;
  
  m = SETWORDSNEEDED(n);
  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

  DYNALLOC1(int, lab, lab_sz, n, "malloc");
  DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
  DYNALLOC1(int, orbits, orbits_sz, n, "malloc");

  densenauty(g, lab, ptn, orbits, &options, &stats, m, n, NULL);
  int group_order;
  group_order = stats.grpsize1;
  int r;
  for (r = 0; r < stats.grpsize2; ++r) group_order *= 10;
    
  group = groupptr(FALSE);

  makecosetreps(group);
  if (target_array_size >= group_order*n)
  {
    store_group(group, target_array);
  }
  else target_array[0] = -2;
}

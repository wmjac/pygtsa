/* Copyright (C) 2014 William M. Jacobs
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "config.h"

#include "sample.h"

histogramEV_t * histogramEV_alloc (const graph_t * graph)
{
  histogramEV_t * h = (histogramEV_t *) malloc (sizeof (histogramEV_t));
  if (h == NULL) goto fail;
  h->max_x = graph->max_nedges + 1;
  h->max_y = graph->max_nedges - graph->max_nvertices + 2;
  h->data = (double *) malloc (sizeof (double) * h->max_x * h->max_y);
  if (h->data == NULL) goto fail;
  return h;
 fail:
  histogramEV_free (h);
  return NULL;
}

void histogramEV_free (histogramEV_t * h)
{
  if (h == NULL) return;
  if (h->data != NULL) free (h->data);
  free (h);
}

void histogramEV_fill (histogramEV_t * h, double z)
{
  size_t i;
  for (i = 0; i < h->max_x * h->max_y; i++)
    {
      h->data[i] = z;
    }
}

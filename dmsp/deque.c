/* movstat/deque.c
 *
 * Double-ended queue based on a circular buffer
 * 
 * Copyright (C) 2018 Patrick Alken
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

/*typedef int deque_type;*/

typedef struct
{
  deque_type *array;
  int head;
  int tail;
  int size;       /* total elements allocated */
} deque;

static deque *deque_alloc(const size_t n);
static void deque_free(deque * d);
static int deque_empty(deque * d);
static int deque_is_empty(const deque * d);
static int deque_is_full(const deque * d);
static int deque_push_front(const deque_type x, deque * d);
static int deque_push_back(const deque_type x, deque * d);
static int deque_pop_front(deque * d);
static int deque_pop_back(deque * d);
static int deque_peek_front(deque_type * x, const deque * d);
static int deque_peek_back(deque_type * x, const deque * d);
static int deque_n(const deque * d);

static deque *
deque_alloc(const size_t n)
{
  deque *d;

  d = calloc(1, sizeof(deque));
  if (d == NULL)
    return NULL;

  d->array = malloc(n * sizeof(deque_type));
  if (d->array == NULL)
    {
      deque_free(d);
      return NULL;
    }

  d->head = -1;
  d->tail = 0;
  d->size = (int) n;

  return d;
}

static void
deque_free(deque * d)
{
  if (d->array)
    free(d->array);

  free(d);
}

/* empty the queue */
static int
deque_empty(deque * d)
{
  d->head = -1;
  d->tail = 0;
  return GSL_SUCCESS;
}

/* check if queue is empty */
static int
deque_is_empty(const deque * d)
{
  return (d->head == -1);
}

/* check if queue is full */
static int
deque_is_full(const deque * d)
{
  return ((d->head == 0 && d->tail == d->size - 1) ||
          (d->head == d->tail + 1));
}

static int
deque_push_front(const deque_type x, deque * d)
{
  if (deque_is_full(d))
    {
      GSL_ERROR("deque is full", GSL_EOVRFLW);
    }
  else
    {
      if (d->head == -1)     /* queue is empty */
        {
          d->head = 0;
          d->tail = 0;
        }
      else if (d->head == 0) /* head is in first position, wrap to end */
        {
          d->head = d->size - 1;
        }
      else                   /* decrement head */
        {
          --(d->head);
        }

      /* insert element */
      d->array[d->head] = x;

      return GSL_SUCCESS;
    }
}

static int
deque_push_back(const deque_type x, deque * d)
{
  if (deque_is_full(d))
    {
      GSL_ERROR("deque is full", GSL_EOVRFLW);
    }
  else
    {
      if (d->head == -1)               /* queue is empty */
        {
          d->head = 0;
          d->tail = 0;
        }
      else if (d->tail == d->size - 1) /* tail is in last position, wrap to 0 */
        {
          d->tail = 0;
        }
      else                             /* increment tail */
        {
          ++(d->tail);
        }

      /* insert element */
      d->array[d->tail] = x;

      return GSL_SUCCESS;
    }
}

static int
deque_pop_front(deque * d)
{
  if (deque_is_empty(d))
    {
      GSL_ERROR("cannot pop element from empty queue", GSL_EOVRFLW);
    }
  else
    {
      if (d->head == d->tail)          /* queue has only one element */
        {
          d->head = -1;
          d->tail = -1;
        }
      else if (d->head == d->size - 1) /* head is in last position, wrap to 0 */
        {
          d->head = 0;
        }
      else                             /* increment head */
        {
          ++(d->head);
        }

      return GSL_SUCCESS;
    }
}

static int
deque_pop_back(deque * d)
{
  if (deque_is_empty(d))
    {
      GSL_ERROR("cannot pop element from empty queue", GSL_EOVRFLW);
    }
  else
    {
      if (d->head == d->tail)    /* queue has only one element */
        {
          d->head = -1;
          d->tail = -1;
        }
      else if (d->tail == 0)     /* tail is in first position, wrap to end */
        {
          d->tail = d->size - 1;
        }
      else                       /* decrement tail */
        {
          --(d->tail);
        }

      return GSL_SUCCESS;
    }
}

static int
deque_peek_front(deque_type * x, const deque * d)
{
  if (deque_is_empty(d))
    {
      GSL_ERROR("queue is empty", GSL_EBADLEN);
    }
  else
    {
      *x = d->array[d->head];
      return GSL_SUCCESS;
    }
}

static int
deque_peek_back(deque_type * x, const deque * d)
{
  if (deque_is_empty(d) || d->tail < 0)
    {
      GSL_ERROR("queue is empty", GSL_EBADLEN);
    }
  else
    {
      *x = d->array[d->tail];
      return GSL_SUCCESS;
    }
}

/* number of items in deque */
static int
deque_n(const deque * d)
{
  if (deque_is_empty(d))
    {
      return 0;
    }
  else
    {
      int size = (d->size - d->head + d->tail) % d->size + 1;
      return size;
    }
}

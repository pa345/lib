/*
 * stack.c
 * Patrick Alken
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

/*typedef double stack_data_t;*/

typedef struct
{
  stack_data_t * array;
  int n;
  int top;
} stack;

static stack *stack_alloc(const int n);
static void stack_free(stack * s);
static int stack_push(const stack_data_t data, stack * s);
static int stack_pop(stack_data_t * data, stack * s);
static int stack_peek(stack_data_t * data, const stack * s);
static int stack_isempty(const stack * s);

/* allocate stack with maximum of n elements of size 'size' */
static stack *
stack_alloc(const int n)
{
  stack *s;

  s = calloc(1, sizeof(stack));
  if (!s)
    {
      fprintf(stderr, "stack_alloc: unable to calloc: %s\n",
              strerror(errno));
      return NULL;
    }

  s->array = malloc(n * sizeof(stack_data_t));
  if (!s->array)
    {
      stack_free(s);
      fprintf(stderr, "stack_alloc: unable to allocate array: %s\n",
              strerror(errno));
      return NULL;
    }

  s->top = -1;
  s->n = n;

  return s;
}

static void
stack_free(stack * s)
{
  if (s->array)
    free(s->array);

  free(s);
}

static int
stack_push(const stack_data_t data, stack * s)
{
  if (s->top >= s->n - 1)
    {
      fprintf(stderr, "stack_push: overflow\n");
      return -1;
    }
  else
    {
      s->array[++(s->top)] = data;
      return 0;
    }
}

static int
stack_pop(stack_data_t * data, stack * s)
{
  if (s->top < 0)
    {
      fprintf(stderr, "stack_pop: underflow\n");
      return -1;
    }
  else
    {
      *data = s->array[(s->top)--];
      return 0;
    }
}

static int
stack_peek(stack_data_t * data, const stack * s)
{
  if (s->top < 0)
    {
      fprintf(stderr, "stack_peek: underflow\n");
      return -1;
    }
  else
    {
      *data = s->array[s->top];
      return 0;
    }
}

static int
stack_isempty(const stack * s)
{
  return (s->top < 0);
}

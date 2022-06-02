/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/network/utilities.h
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <stdlib.h>

#include "arepoconfig.h"

#if !defined(__has_attribute) && !defined(__GNUC__)
/* remove GCC-style attributes if the compiler does not support them */
#define __attribute__(x)
#endif

#define safe_fgets(str, num, stream) util_fgets(str, num, stream, (char *)__FILE__, __LINE__)
#define safe_fread(ptr, size, count, stream) util_fread(ptr, size, count, stream, (char *)__FILE__, __LINE__)

void myprintf(const char *format, ...) __attribute__((format(printf, 1, 2)));
char *util_fgets(char *str, int num, FILE *stream, char *file, int line);
size_t util_fread(void *ptr, size_t size, size_t count, FILE *stream, char *file, int line);

#endif /* UTILITIES_H */

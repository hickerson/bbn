/* testF.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int fortfunc_(integer *ii, real *ff)
{
    /* Format strings */
    static char fmt_100[] = "(\002ii=\002,i2,\002 ff=\002,f6.3)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, fmt_100, 0 };


    s_wsfe(&io___1);
    do_fio(&c__1, (char *)&(*ii), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*ff), (ftnlen)sizeof(real));
    e_wsfe();
    return 0;
} /* fortfunc_ */


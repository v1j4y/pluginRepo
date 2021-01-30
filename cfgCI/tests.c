int main() {
#include <stdio.h>
#include "munit/munit.h"

void tests(){

    int foo = 1;
    int bar = 1;

    munit_assert_int(foo, ==, bar);

    printf("Hello world (munit-testing)\n");
}

return 0;
}

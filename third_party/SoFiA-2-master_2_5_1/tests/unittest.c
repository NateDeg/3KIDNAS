#include <stdlib.h>
#include "test_LinkerPar.h"

// Run unittest suite
int main(void) {
    int no_failed = 0;                   
    Suite *s;
    SRunner *runner;                     

    s = suite_create("SoFiA-2");
    runner = srunner_create(s);
    srunner_add_suite(runner, LinkerPar_test_suite());

    srunner_run_all(runner, CK_NORMAL);  
    no_failed = srunner_ntests_failed(runner); 
    srunner_free(runner);                      
    return (no_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;  
}
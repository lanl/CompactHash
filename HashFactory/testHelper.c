/* Copyright 2013-14.  Los Alamos National Security, LLC. This material was produced
 * under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
 * Laboratory (LANL), which is operated by Los Alamos National Security, LLC
 * for the U.S. Department of Energy. The U.S. Government has rights to use,
 * reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS
 * ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified
 * to produce derivative works, such modified software should be clearly marked,
 * so as not to confuse it with the version available from LANL.   
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy
 * of the License at 
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed
 * under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 *
 * Under this license, it is required to include a reference to this work. We
 * request that each derivative work contain a reference to LANL Copyright 
 * Disclosure C14043/LA-CC-14-003 so that this work's impact can be roughly
 * measured. In addition, it is requested that a modifier is included as in
 * the following example:
 *
 * //<Uses | improves on | modified from> LANL Copyright Disclosure C14043/LA-CC-14-003
 *
 * This is LANL Copyright Disclosure C14043/LA-CC-14-003
 */

/**
 * @file   test.c
 * @author Peter Ahrens
 * @date   Fri Jun 7 2013 
 */

#include "testHelper.h"

void Test_StartTry(char *nameOfTest){
  printf("TRYING %s: ", nameOfTest);
  fflush(stdout);
}
    
void Test_EndTry(TestData *currentTest){
  printf("PASSED\n");
  currentTest->numPassed++;
}

void Test_PointerEquality(TestData *currentTest, char *nameOfTest, void *expected, void *observed, char *failMessage){
  if(expected == observed){
    printf("%s: PASSED\n", nameOfTest);
    currentTest->numPassed++;
  }else{
    printf("%s: FAILED\n", nameOfTest);
    printf("  POINTERS DO NOT MATCH");
    if(failMessage != NULL){
      printf("  ADDITIONAL INFO: %s\n", failMessage);
    }
    currentTest->numFailed++;
    printf("\n");
  }
}

void Test_PointerInequality(TestData *currentTest, char *nameOfTest, void *notExpected, void *observed, char *failMessage){
  if(notExpected != observed){
    printf("%s: PASSED\n", nameOfTest);
    currentTest->numPassed++;
  }else{
    printf("%s: FAILED\n", nameOfTest);
    printf("  POINTERS SHOULD NOT MATCH");
    if(failMessage != NULL){
      printf("  ADDITIONAL INFO: %s\n", failMessage);
    }
    currentTest->numFailed++;
    printf("\n");
  }
}

void Test_IntEquality(TestData *currentTest, char *nameOfTest, int expected, int observed, char *failMessage){
  if(expected == observed){
    printf("%s: PASSED\n", nameOfTest);
    currentTest->numPassed++;
  }else{
    printf("%s: FAILED\n", nameOfTest);
    printf("  EXPECTED: %d\n", expected);
    printf("  OBSERVED: %d\n", observed);
    if(failMessage != NULL){
      printf("  ADDITIONAL INFO: %s\n", failMessage);
    }
    currentTest->numFailed++;
    printf("\n");
  }
}

void Test_IntArrayEquality(TestData *currentTest, char *nameOfTest, unsigned int size, int *expected, int *observed, char *failMessage){
  int result = 0;
  int i;
  for(i = 0; i < size; i++){
    if(expected[i] != observed[i]){
      result = 1;
    }
  }
  if(result == 0){
    printf("%s: PASSED\n", nameOfTest);
    currentTest->numPassed++;
  }else{
    printf("%s: FAILED\n", nameOfTest);
    printf("  EXPECTED: {");
    for(i = 0; i < size; i++){
      printf("%d, ", expected[i]);
    }
    printf("}\n");
    printf("  OBSERVED: {");
    for(i = 0; i < size; i++){
      printf("%d, ", observed[i]);
    }
    printf("}\n");
    if(failMessage != NULL){
      printf("  ADDITIONAL INFO: %s\n", failMessage);
    }
    currentTest->numFailed++;
    printf("\n");
  }
}

void Test_IntLessThan(TestData *currentTest, char *nameOfTest, int expectedAbove, int observed, char *failMessage){
  if(expectedAbove > observed){
    printf("%s: PASSED\n", nameOfTest);
    currentTest->numPassed++;
  }else{
    printf("%s: FAILED\n", nameOfTest);
    printf("  EXPECTED OBSERVED VALUE TO BE LESS THAN: %d\n", expectedAbove);
    printf("  OBSERVED: %d\n", observed);
    if(failMessage != NULL){
      printf("  ADDITIONAL INFO: %s\n", failMessage);
    }
    currentTest->numFailed++;
    printf("\n");
  }
}

void Test_IntGreaterThan(TestData *currentTest, char *nameOfTest, int expectedBelow, int observed, char *failMessage){
  if(expectedBelow < observed){
    printf("%s: PASSED\n", nameOfTest);
    currentTest->numPassed++;
  }else{
    printf("%s: FAILED\n", nameOfTest);
    printf("  EXPECTED OBSERVED VALUE TO BE GREATER THAN: %d\n", expectedBelow);
    printf("  OBSERVED: %d\n", observed);
    if(failMessage != NULL){
      printf("  ADDITIONAL INFO: %s\n", failMessage);
    }
    currentTest->numFailed++;
    printf("\n");
  }
}

void Test_FinishSubTest(TestData *currentTest, TestData *subTest){
  if(subTest->numFailed == 0){
    currentTest->numPassed++;
  }else{
    currentTest->numFailed++;
  }
  Test_FinishTest(subTest);
}
  

TestData *Test_CreateTest(char *nameOfTestGroup){
  TestData *newTest = (TestData*)malloc(sizeof(TestData));
  newTest->nameOfTestGroup = nameOfTestGroup;
  newTest->numPassed = 0;
  newTest->numFailed = 0;
  printf("%s STARTED\n", nameOfTestGroup);
  return newTest;
}

void Test_FinishTest(TestData *currentTest){
  printf("%s FINISHED\n", currentTest->nameOfTestGroup);
  printf("%d%% (%d/%d) TESTS PASSED\n\n", 
         ((currentTest->numPassed * 100) / (currentTest->numPassed + currentTest->numFailed)),
         currentTest->numPassed,
         currentTest->numPassed + currentTest->numFailed);
  free(currentTest);
}

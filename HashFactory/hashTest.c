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
 * @file   hashTest.c
 * @author Peter Ahrens
 * @date   Fri Jun 7 2013 
 */

#include "testHelper.h"
#include "HashFactory.h"

void sanityTest(int hashType, TestData *metaTest, char* subTestName){
  TestData *sanityTestData = Test_CreateTest(subTestName);
  int emptyValue = -1;
  intintHash_Factory *sanityFactory = intintHash_CreateFactory(hashType, &emptyValue, 128, NULL, NULL);
  intintHash_Table *sanityTable = intintHash_CreateTable(sanityFactory, hashType, 21, 7, 0.5);
  Test_IntEquality(sanityTestData, "Correct Hash Type Check", hashType, intintHash_GetTableType(sanityTable), NULL);
  intintHash_EmptyTable(sanityTable);
  int keys[] = {6, 6, 9, 12, 15, 18, 21};
  int values[] = {101, 101, 102, 103, 104, 105, 106};
  int keys1[] = {6, 9, 12, 15, 18, 21};
  int values1[] = {101, 102, 103, 104, 105, 106};
  int values2[] = {201, 202, 203, 204, 205, 206};
  int valuesOut[6];
  intintHash_Insert(sanityTable, 7, keys, values);
  int sanityInt;
  intintHash_QuerySingle(sanityTable, 6, &sanityInt);
  Test_IntEquality(sanityTestData, "QuerySingle after Insert Check 1", 101, sanityInt, NULL);
  intintHash_QuerySingle(sanityTable, 9, &sanityInt);
  Test_IntEquality(sanityTestData, "QuerySingle after Insert Check 2", 102, sanityInt, NULL);
  intintHash_QuerySingle(sanityTable, 12, &sanityInt);
  Test_IntEquality(sanityTestData, "QuerySingle after Insert Check 3", 103, sanityInt, NULL);
  intintHash_QuerySingle(sanityTable, 15, &sanityInt);
  Test_IntEquality(sanityTestData, "QuerySingle after Insert Check 4", 104, sanityInt, NULL);
  intintHash_QuerySingle(sanityTable, 18, &sanityInt);
  Test_IntEquality(sanityTestData, "QuerySingle after Insert Check 5", 105, sanityInt, NULL);
  intintHash_QuerySingle(sanityTable, 21, &sanityInt);
  Test_IntEquality(sanityTestData, "QuerySingle after Insert Check 6", 106, sanityInt, NULL);
  intintHash_InsertSingleNoOverwrite(sanityTable, 6, 201);
  intintHash_InsertSingleNoOverwrite(sanityTable, 9, 202);
  intintHash_InsertSingleNoOverwrite(sanityTable, 12, 203);
  intintHash_InsertSingleNoOverwrite(sanityTable, 15, 204);
  intintHash_InsertSingleNoOverwrite(sanityTable, 18, 205);
  intintHash_InsertSingleNoOverwrite(sanityTable, 21, 206);
  intintHash_Query(sanityTable, 6, keys1, &valuesOut[0]);
  Test_IntArrayEquality(sanityTestData, "Query after InsertSingleNoOverwrite Check", 6, values1, valuesOut, NULL);
  intintHash_InsertSingle(sanityTable, 6, 201);
  intintHash_InsertSingle(sanityTable, 9, 202);
  intintHash_InsertSingle(sanityTable, 12, 203);
  intintHash_InsertSingle(sanityTable, 15, 204);
  intintHash_InsertSingle(sanityTable, 18, 205);
  intintHash_InsertSingle(sanityTable, 21, 206);
  intintHash_Query(sanityTable, 6, keys1, valuesOut);
  Test_IntArrayEquality(sanityTestData, "Query after InsertSingle Check", 6, values2, valuesOut, NULL);
  intintHash_InsertNoOverwrite(sanityTable, 6, keys1, values1);
  intintHash_Query(sanityTable, 6, keys1, valuesOut);
  Test_IntArrayEquality(sanityTestData, "Query after InsertNoOverwrite Check", 6, values2, valuesOut, NULL);
  intintHash_DestroyTable(sanityTable);
  intintHash_DestroyFactory(sanityFactory);
  Test_FinishSubTest(metaTest, sanityTestData);
}

int main(int argc, char *argv[]){
  CLHash_Init(argv[0]);
  TestData *allHashTestsData = Test_CreateTest("HashTest");
  sanityTest(IDENTITY_PERFECT_HASH_ID, allHashTestsData, "intintIdentityPerfectHash Test");
  sanityTest(IDENTITY_SENTINEL_PERFECT_HASH_ID, allHashTestsData, "intintIdentitySentinelPerfectHash Test");
  sanityTest(LCG_LINEAR_OPEN_COMPACT_HASH_ID, allHashTestsData, "intintLCGLinearOpenCompactHash Test");
  sanityTest(LCG_QUADRATIC_OPEN_COMPACT_HASH_ID, allHashTestsData, "intintLCGQuadraticOpenCompactHash Test");
  sanityTest(IDENTITY_PERFECT_CL_HASH_ID, allHashTestsData, "intintIdentityPerfectCLHash Test");
  sanityTest(IDENTITY_SENTINEL_PERFECT_CL_HASH_ID, allHashTestsData, "intintIdentitySentinelPerfectCLHash Test");
  sanityTest(LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID, allHashTestsData, "intintLCGLinearOpenCompactCLHash Test");
  sanityTest(LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID, allHashTestsData, "intintLCGQuadraticOpenCompactCLHash Test");
  Test_FinishTest(allHashTestsData);
  return 0;
}

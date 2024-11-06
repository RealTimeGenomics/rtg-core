/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.variant.bayes.multisample;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.util.ArrayList;

import com.rtg.util.PortableRandom;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.Variant;

import junit.framework.TestCase;

/**
 * Test class
 */
public class BedComplexitiesWriterTest extends TestCase {
  private byte[] template(int length) {
    final PortableRandom r = new PortableRandom(42);
    final byte[] res = new byte[length];
    for (int i = 0; i < length; ++i) {
      res[i] = (byte) r.nextInt(5);
    }
    return res;
  }

  public void test() throws IOException {
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(6));
    chunk.add(TestUtils.createVariant(10));
    chunk.add(TestUtils.createVariant(15));

    final Complexities regions = new Complexities(chunk, "foo", 0, 50, 3, 15, template(30), true, null);
    Complexities.fixDangling(null, regions);
    Complexities.fixDangling(regions, null);
    final MemoryPrintStream ms = new MemoryPrintStream();
    final BedComplexitiesWriter bed = new BedComplexitiesWriter(ms.outputStream(), "foobar", 0);
    bed.write(regions);
    bed.finish();
    assertEquals("foobar\t4\t7\tcomplex-called" + LS
        + "foobar\t10\t16\tcomplex-called" + LS, ms.toString());
  }

  public void testMixture() throws IOException {
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(6));
    chunk.add(TestUtils.createVariant(10));
    chunk.add(TestUtils.createVariant(12));
    chunk.add(TestUtils.createVariant(14));
    chunk.add(TestUtils.createVariant(15));
    chunk.add(TestUtils.createVariant(17));
    chunk.add(TestUtils.createVariant(19));
    chunk.add(TestUtils.createVariant(30, true));

    final Complexities regions = new Complexities(chunk, "foo", 0, 50, 3, 7, template(34), true, null);
    Complexities.fixDangling(null, regions);
    Complexities.fixDangling(regions, null);
    final MemoryPrintStream ms = new MemoryPrintStream();
    final BedComplexitiesWriter bed = new BedComplexitiesWriter(ms.outputStream(), "foobar", 0);
    bed.write(regions);
    bed.finish();
    assertEquals(""
        + "foobar\t4\t7\tcomplex-called" + LS
        + "foobar\t10\t20\thyper-complex" + LS
        + "foobar\t30\t31\tcomplex-called" + LS, ms.toString());
  }

  // test multiple complexity spanning hyper complex region
  public void testHyperComplexMultispan() throws IOException {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(13));
    chunkA.add(TestUtils.createVariant(15));
    chunkA.add(TestUtils.createVariant(17));
    chunkA.add(TestUtils.createVariant(19));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 5, template(30), true, null);
    final ArrayList<Variant> chunkB = new ArrayList<>();
    for (int i = 21; i < 40; i += 2) {
      chunkB.add(TestUtils.createVariant(i, true));
    }
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 5, template(30), true, null);

    final ArrayList<Variant> chunkC = new ArrayList<>();
    chunkC.add(TestUtils.createVariant(41));
    chunkC.add(TestUtils.createVariant(43));
    chunkC.add(TestUtils.createVariant(45));
    chunkC.add(TestUtils.createVariant(47));
    final Complexities regionsC = new Complexities(chunkC, "foo", 40, 60, 3, 5, template(30), true, null);
    regionsC.globalIntegrity();
    //merge a and b
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    //merge b and c
    Complexities.fixDangling(regionsB, regionsC);
    Complexities.fixDangling(regionsC, null);
    final MemoryPrintStream ms = new MemoryPrintStream();
    final BedComplexitiesWriter bed = new BedComplexitiesWriter(ms.outputStream(), "foobar", 0);
    bed.write(regionsA);
    bed.write(regionsB);
    bed.write(regionsC);
    bed.finish();
    assertEquals("foobar\t13\t48\thyper-complex" + LS, ms.toString());
  }

  // test hyper complex spanning single boundary
  public void testHyperComplexSingleSpan() throws IOException {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(13));
    chunkA.add(TestUtils.createVariant(15));
    chunkA.add(TestUtils.createVariant(17));
    chunkA.add(TestUtils.createVariant(19));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 5, template(30), true, null);
    final ArrayList<Variant> chunkB = new ArrayList<>();
    for (int i = 21; i < 30; i += 2) {
      chunkB.add(TestUtils.createVariant(i, true));
    }
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 5, template(30), true, null);

    //merge a and b
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    final MemoryPrintStream ms = new MemoryPrintStream();
    final BedComplexitiesWriter bed = new BedComplexitiesWriter(ms.outputStream(), "foobar", 0);
    bed.write(regionsA);
    bed.write(regionsB);
    bed.finish();
    assertEquals("foobar\t13\t30\thyper-complex" + LS, ms.toString());
  }

  // test hyper complex spanning single boundary twice
  public void testTwoHyperComplexSingleSpan() throws IOException {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(13));
    chunkA.add(TestUtils.createVariant(15));
    chunkA.add(TestUtils.createVariant(17));
    chunkA.add(TestUtils.createVariant(19));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 5, template(30), true, null);
    final ArrayList<Variant> chunkB = new ArrayList<>();
    for (int i = 21; i < 30; i += 2) {
      chunkB.add(TestUtils.createVariant(i, true));
    }
    chunkB.add(TestUtils.createVariant(39));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 5, template(30), true, null);

    final ArrayList<Variant> chunkC = new ArrayList<>();
    chunkC.add(TestUtils.createVariant(41));
    chunkC.add(TestUtils.createVariant(43));
    chunkC.add(TestUtils.createVariant(45));
    final Complexities regionsC = new Complexities(chunkC, "foo", 40, 60, 3, 5, template(30), true, null);

    //merge a and b
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, regionsC);
    Complexities.fixDangling(regionsC, null);
    final MemoryPrintStream ms = new MemoryPrintStream();
    final BedComplexitiesWriter bed = new BedComplexitiesWriter(ms.outputStream(), "foobar", 0);
    bed.write(regionsA);
    bed.write(regionsB);
    bed.write(regionsC);
    bed.finish();
    assertEquals(""
        + "foobar\t13\t30\thyper-complex" + LS
        + "foobar\t39\t46\thyper-complex" + LS, ms.toString());
  }

  //simple adjacent interesting calls on boundary - there is a complex object dupicated in the two complexity objects the second should be ignored
  public void testSpanBoundary() throws IOException {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(19));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), true, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(21));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 15, template(30), true, null);

    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);

    final MemoryPrintStream ms = new MemoryPrintStream();
    final BedComplexitiesWriter bed = new BedComplexitiesWriter(ms.outputStream(), "foobar", 0);
    bed.write(regionsA);
    bed.write(regionsB);
    bed.finish();
    assertEquals("foobar\t19\t22\tcomplex-called" + LS, ms.toString());
  }

  public void testIndelOnBoundary() throws IOException {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(20, 20, true, 0));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), true, null);
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, null);
    final MemoryPrintStream ms = new MemoryPrintStream();
    final BedComplexitiesWriter bed = new BedComplexitiesWriter(ms.outputStream(), "foobar", 0);
    bed.write(regionsA);
    bed.finish();
    assertEquals("foobar\t20\t20\tcomplex-called" + LS, ms.toString());
  }

  public void testOverflowBoundary() throws IOException {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    final Variant varA = TestUtils.createVariant(14896, 14999, true, 0);
    varA.setOverflow();
    chunkA.add(varA);
    final Variant varB = TestUtils.createVariant(15000, 15001, false, 0);
    varB.setInteresting();
    varB.setNonIdentityPosterior(0.612951146759638);
    chunkB.add(varB);
    final Complexities regionsA = new Complexities(chunkA, "GL000226.1", 0, 15000, 3, 15, template(15008), false, null);
    final Complexities regionsB = new Complexities(chunkB, "GL000226.1", 15000, 15008, 3, 15, template(15008), false, null);
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    final MemoryPrintStream ms = new MemoryPrintStream();
    final BedComplexitiesWriter bed = new BedComplexitiesWriter(ms.outputStream(), "GL000226.1", 0);
    bed.write(regionsA);
    bed.write(regionsB);
    bed.finish();
    assertEquals("GL000226.1\t14896\t14999\textreme-coverage" + LS, ms.toString());
  }
}

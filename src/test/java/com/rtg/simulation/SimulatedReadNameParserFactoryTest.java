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
package com.rtg.simulation;

import junit.framework.TestCase;

/**
 * Test class
 */
public class SimulatedReadNameParserFactoryTest extends TestCase {

  public void testFactory() {
    SimulatedReadNameParser p = SimulatedReadNameParserFactory.getParser("read9 0/0/simulatedSequence1/9068/F/7.1X27.");
    assertNotNull(p);
    assertTrue(p.getClass().equals(NewReadNameParser.class));

    p = SimulatedReadNameParserFactory.getParser("read9/0/0/simulatedSequence1/9068/F/7.1X27.");
    assertNotNull(p);
    assertTrue(p.getClass().equals(NewReadNameParser.class));

    p = SimulatedReadNameParserFactory.getParser("read02:simulatedSequence1:30713R:S0:I0:D0");
    assertNotNull(p);
    assertTrue(p instanceof OldReadNameParser);

    p = SimulatedReadNameParserFactory.getParser("HISEQ1:18:H8VC6ADXX:2:2104:14516:42086");
    assertNull(p);

    p = SimulatedReadNameParserFactory.getParser("0 frag9/0/0/simulatedSequence1/9068/F/7.1X27.");
    assertNotNull(p);
    assertTrue(p instanceof NewestReadNameParser);

    p = SimulatedReadNameParserFactory.getParser("@CFTR.3.70s_282_148_0_1_0_0_4:0:0_4:0:0_10c/1");
    assertNotNull(p);
    assertTrue(p instanceof DwgsimReadNameParser);

    assertNull(SimulatedReadNameParserFactory.getParser("cgread2:si6"));
    assertNull(SimulatedReadNameParserFactory.getParser("cgread2/si6"));
    assertNull(SimulatedReadNameParserFactory.getParser("/////"));
    assertNull(SimulatedReadNameParserFactory.getParser("::::"));
    assertNull(SimulatedReadNameParserFactory.getParser(""));
  }
}

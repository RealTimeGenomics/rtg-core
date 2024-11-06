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
package com.rtg.ml;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class MlDataTypeTest extends TestCase {
  public void test() throws IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    try (final DataOutputStream dos = new DataOutputStream(mps.outputStream())) {
      MlDataType.DOUBLE.save(5.5, dos);
      MlDataType.BOOLEAN.save(true, dos);
      MlDataType.INTEGER.save(556, dos);
      MlDataType.STRING.save("foorahgah", dos);
    }
    try (final DataInputStream dis = new DataInputStream(new ByteArrayInputStream(mps.toByteArray()))) {
      assertEquals(5.5, MlDataType.DOUBLE.load(dis));
      assertEquals(true, MlDataType.BOOLEAN.load(dis));
      assertEquals(556, MlDataType.INTEGER.load(dis));
      assertEquals("foorahgah", MlDataType.STRING.load(dis));
    }
  }
}

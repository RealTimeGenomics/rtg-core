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

package com.rtg.assembler.graph.io;

import java.io.IOException;
import java.text.ParseException;
import java.util.Collections;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.UUID;

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.implementation.GraphImplementation;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.Talkback;
import com.rtg.util.store.StoreDirString;
import com.rtg.util.store.StoreDirectory;
import com.rtg.util.store.StoreFile;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 */
public class GraphWriterTest extends TestCase {

  protected NanoRegression mNano = null;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass(), false);
  }

  @Override
  public void tearDown() throws Exception {
    // clear the module name so later tests don't report SlimException to the
    // Talkback system
    Talkback.setModuleName(null);
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  public void testHeader() throws ParseException, IOException {
    final Map<String, String> contigAttributes = new LinkedHashMap<>();
    contigAttributes.put("con1", "comment1");
    contigAttributes.put("con2", "comment2");
    final Map<String, String> pathAttributes = new LinkedHashMap<>();
    pathAttributes.put("path1", "pcomment1");
    pathAttributes.put("path2", "pcomment2");
    final Set<UUID> inputUids = new LinkedHashSet<>();
    inputUids.add(new UUID(21, 103));
    inputUids.add(new UUID(1, 1));

    final Graph graph = new GraphImplementation(0, contigAttributes, pathAttributes);
    final Date date = GraphWriter.dateFormat().parse("2012/05/15 16:03:00");
    final UUID uuid = new UUID(42, 101);
    final StoreDirectory dir = new StoreDirString();
    final GraphWriter gw = new GraphWriter(graph, dir, 0, false);
    gw.writeHeader(graph, "test", date, uuid, inputUids);
    final String actual = dir.child("header.tsv").content();
    //System.err.println(actual);
    mNano.check("header2", actual);
  }

  public void testLongContig() throws IOException {
    final GraphImplementation graph = new GraphImplementation(0, Collections.emptyMap(), Collections.emptyMap());
    final StringBuilder sb = new StringBuilder();
    final String tenLong = "ACTG" + "GGTT" + "AA";
    for (int i = 0; i < 16; ++i) {
      sb.append(tenLong);
    }

    final String contigStr = sb.toString();
    assertEquals(160, contigStr.length());
    final Contig contig = new ContigString(contigStr);
    graph.addContig(contig);
    graph.addContig(new ContigString("AAAA"));
    graph.addContig(new ContigString("CCCC"));

    final StoreDirectory dir = new StoreDirString();
    GraphWriter.write(graph, dir, "test", Collections.emptySet());

    final StoreFile contigs = dir.child("contig.1.fa");
    //System.err.print(contigs.content());
    //System.err.println("====");
    mNano.check("contigLong", contigs.content());

    final StoreFile paths = dir.child("path.1.tsv");
    //System.err.print(paths.content());
    mNano.check("pathsLong", paths.content());
    //System.err.println("====");

    final StoreFile header = dir.child("header.tsv");
    //System.err.print(header.content());
    //System.err.println("====");

    final String rawContent = header.content();
    final String editedContent = StringUtils.grepMinusV(StringUtils.grepMinusV(rawContent, "date"), "guid");
    //System.err.println(editedContent);
    mNano.check("headerLong", editedContent);
  }
}

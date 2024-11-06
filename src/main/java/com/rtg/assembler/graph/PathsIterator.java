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

package com.rtg.assembler.graph;

/**
 * Holds the result of a search for all the paths associated with a contig in a <b>Graph</b>
 * (see <code>PathIterator paths(long contigId, PathIterator iter)</code>).
 * <p>
 * Example code using this class:
 * <pre>
 * PathIterator iter = ....;
 * ......
 * iter = graph.paths(contigId, iter);
 * while (true) {
 *    final long path = iter.nextPathId();
 *    if (pathId == 0) {
 *      break;
 *    }
 *    final int length = graph.pathLength(path);
 *    final int index = iter.contigIndex();
 *    assert 0 &lt;= index &amp;&amp; index &lt; length;
 *    assert graph.pathContig(path, index) == contigId;
 * }
 * </pre>
 */
public interface PathsIterator {

  /**
   * @return the graph this iterator operates over.
   */
  Graph graph();

  /**
   * Reset the iterator to all paths that include the <code>contigId</code> (the direction of the path is
   * determined by the sign of <code>contigId</code> ).
   * Equivalent to <code>set(contigId, false)</code>.
   * @param contigId directional contig identifier.
   * @throws IllegalArgumentException if <code>abs(contigId)</code> is 0 or &gt; <code>graph().numberContigs()</code>.
   */
  void set(long contigId);

  /**
   * Reset the iterator to all paths that include the <code>contigId</code> (the direction of the path is
   * determined by the sign of <code>contigId</code> ).
   * If <code>showDeleted</code> is true then all deleted paths are included in the result, otherwise only current (non-deleted) paths are included.
   * @param contigId directional contig identifier.
   * @param showDeleted if true then include deleted paths.
   * @throws IllegalArgumentException if <code>abs(contigId)</code> is 0 or &gt; <code>graph().numberContigs()</code>.
   */
  void set(long contigId, boolean showDeleted);

  /**
   * Get the index within the current path (as given by the last call to <code>nextPathId()</code>) of the contig used to select the paths.
   *
   * @return the index  within the current path of the contig used to set the <code>PathIterator</code>.
   * @throws IllegalStateException if <code>nextPathId()</code> has not been called since the last
   */
  int contigIndex();

  /**
   * Once this has returned a 0 subsequent calls will return 0 until it is reset.
   * @return the path identifier of the next path. 0 if there are no more paths.
   */
  long nextPathId();
}

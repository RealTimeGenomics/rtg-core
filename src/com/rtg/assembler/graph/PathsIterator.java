/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
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
 *    assert 0 <= index && index < length;
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
  int contigIndex() throws IllegalStateException;


  /**
   * Once this has returned a 0 subsequent calls will return 0 until it is reset.
   * @return the path identifier of the next path. 0 if there are no more paths.
   */
  long nextPathId();
}

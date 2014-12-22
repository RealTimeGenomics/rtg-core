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

import java.util.Map;

/**
 * Graph which has the notion of contigs (a nucleotide sequence with a direction), paths through the contigs, and
 * the ability to mark nodes as deleted but optionally retained in the graph.
 * <p>
 * Representations of various atomic values:
 * <ul>
 * <li>nucleotide (nt) - N=0, A=1, C=2, G=3, T=4 </li>
 * <li><code>contigId</code> - long &ne; 0, a contig and its direction. 0 reserved for a null non-existent contig.
 * The absolute value uniquely represents the contig and the sign the direction (positive is forward through the contig).
 * Allocated in sequence from 1 up.</li>
 * <li><code>pathId</code> - long &ne; 0, a path and its direction. 0 reserved for a null non-existent path.
 * The absolute value uniquely represents the path and the sign the direction.
 * Allocated in sequence from 1 up.</li>
 * </ul>
 */
public interface Graph {

  /**
   * Get the number of contigs. All <code>abs(contigId)</code>s are &le; this value.
   * @return the number of contigs in the graph including ones that are marked as deleted.
   */
  long numberContigs();

  /**
   * Get the number of paths. All <code>abs(pathId)</code>s are &le; this value.
   * @return the number of contigs in the graph including ones that are marked as deleted.
   */
  long numberPaths();

  /**
   * Check if a contig has been deleted.
   * @param contigId identifier of contig being checked.
   * @return true iff this contig has been deleted.
   * @throws IllegalArgumentException if <code>abs(contigId)</code> is 0 or &gt; <code>numberContigs()</code>.
   */
  boolean contigDeleted(long contigId) throws IllegalArgumentException;

  //Nucleotides

  /**
   * Get the length of the contig.
   * @param contigId directional identifier of contig being checked.
   * @return the number of nucleotides in the contig.
   * @throws IllegalArgumentException if <code>abs(contigId)</code> is 0 or &gt; <code>numberContigs()</code>.
   */
  int contigLength(long contigId) throws IllegalArgumentException;


  /**
   * Get a nucleotide from inside the contig. If the <code>contigId</code> is negative then the nucleotide returned
   * if the reverse complement of the positive version and counts from the other end.
   * @param contigId directional identifier of contig.
   * @param index position selecting the nucleotide (<code>0 &le; index &lt; lengthContig(contigId)</code>).
   * @return the <code>index'th</code> nucleotide in the contig.
   * @throws IllegalArgumentException if <code>abs(contigId)</code> is 0 or &gt; <code>numberContigs()</code>.
   * @throws RuntimeException if <code>index &lt; 0</code> or <code>index &ge; lengthContig(contigId)</code>.
   */
  byte nt(long contigId, int index) throws IllegalArgumentException;

  /**
   * Extract contig information into a <code>Contig</code>.
   * This a slow way (but possibly more object oriented and clearer) of accessing information about a contig.
   * Use <code>lengthContig(contigId)</code> and <code>nt(contigId, index)</code> for fast access.
   * @param contigId directional contig identifier.
   * @return a <code>Contig</code> that holds the information associated with <code>contigId</code>.
   * @throws IllegalArgumentException if <code>abs(contigId)</code> is 0 or &gt; <code>numberContigs()</code>.
   */
  Contig contig(long contigId) throws IllegalArgumentException;

  //Attributes

  /**
   * @return the keys and comments for all the optional attributes available for contigs.
   */
  Map<String, String> contigAttributes();

  /**
   * Get an attribute for a contig. The result is independent of the sign (direction) of <code>contigId</code>.
   * @param contigId directional identifier of contig.
   * @param attribute key for an attribute.
   * @return the attribute value (null if not stored or there is no such attribute).
   * @throws IllegalArgumentException if <code>abs(contigId)</code> is 0 or &gt; <code>numberContigs()</code>.
   */
  String contigAttribute(long contigId, String attribute);

  /**
   * @return the keys and comments for all the optional attributes available for paths.
   */
  Map<String, String> pathAttributes();

  /**
   * Get an attribute for a path. The result is independent of the sign (direction) of <code>pathId</code>.
   * @param pathId identifier of path.
   * @param attribute key for an attribute.
   * @return the attribute value (null if not stored or there is no such attribute).
   * @throws IllegalArgumentException if <code>pathId</code> is 0.
   */
  String pathAttribute(long pathId, String attribute) throws IllegalArgumentException;

  /**
   * Check if a path has been deleted.
   * @param pathId directional identifier of path.
   * @return true iff this path has been deleted.
   * @throws IllegalArgumentException if <code>abs(pathId)</code> is 0 or &gt; <code>numberPaths()</code>.
   */
  boolean pathDeleted(long pathId) throws IllegalArgumentException;

  /**
   * This enables a faster way of finding paths that does not require reconstructing <code>PathIterator</code>s for each call.
   * @return an iterator over this graph for use in the <code>paths(...)</code> calls.
   */
  PathsIterator iterator();

  /**
   * Construct a new iterator over all current (non-deleted) paths that include the <code>contigId</code> (the direction of the path is
   * determined by the sign of <code>contigId</code>).
   * Equivalent to <code>paths(contigId, false)</code>.
   * @param contigId directional contig identifier.
   * @return iterator over the the paths.
   * @throws IllegalArgumentException if <code>abs(contigId)</code> is 0 or &gt; <code>numberContigs()</code>.
   */
  PathsIterator paths(long contigId) throws IllegalArgumentException;


  /**
   * Construct a new iterator over all current (non-deleted) paths that include the <code>contigId</code> (the direction of the path is
   * determined by the sign of <code>contigId</code> ).
   * If <code>showDeleted</code> is true then all deleted paths are included in the result, otherwise only current (non-deleted) paths are included.
   * @param contigId directional contig identifier.
   * @param showDeleted if true then include deleted paths.
   * @return iterator over all the selected paths.
   * @throws IllegalArgumentException if <code>abs(contigId)</code> is 0 or &gt; <code>numberContigs()</code>.
   */
  PathsIterator paths(long contigId, boolean showDeleted) throws IllegalArgumentException;

  /**
   * Get the number of contigs included in a path.
   * @param pathId directional path identifier.
   * @return the number of contigs (&ge; 2).
   * @throws IllegalArgumentException if <code>pathId</code> is 0.
   */
  int pathLength(long pathId) throws IllegalArgumentException;

  /**
   * Get one contig from the path.
   * @param pathId directional identifier of path.
   * @param index select the contig on the path (0 based, <code>0 &le; index &lt; pathLength(pathId)</code>).
   * @return a directional contig identifier (the direction is determined by the overall direction of the path and the direction of the contig within the path).
   * @throws IllegalArgumentException if <code>abs(pathId)</code> is 0 or &gt; <code>numberPaths()</code>.
   * @throws RuntimeException <code>index &lt; 0</code> or <code>index &ge; pathLength(pathId)</code>.
   */
  long pathContig(long pathId, int index) throws IllegalArgumentException;

  /**
   * Extract path information into a <code>Path</code>.
   * This a slow way (but possibly more object oriented and clearer) of accessing information about a path.
   * Use <code>pathLength(pathId)</code> and <code>pathContig(pathId, index)</code> for fast access.
   * @param pathId directional identifier of path.
   * @return a <code>Path</code> that holds the information associated with <code>pathId</code>.
   * @throws IllegalArgumentException if <code>abs(pathId)</code> is 0 or &gt; <code>numberPaths()</code>.
   */
  Path path(long pathId) throws IllegalArgumentException;

  /**
   * @return the number of bases adjacent contigs overlap within the graph
   */
  int contigOverlap();
}

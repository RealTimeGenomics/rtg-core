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
 */
public interface MutableGraph extends Graph {

  /**
   * Add another contig to the graph copying its contents from contig and allocating and
   * returning a new unique contig identifier.
   * @param contig the contig used to set the nucleotides stored.
   * @return the unique identifier for the contig (&gt; 0).
   * @throws RuntimeException if <code>contig</code> is null.
   */
  long addContig(Contig contig);

  /**
   * Mark this contig as as deleted. However, see the <code>paths(int contigId, PathIterator iter, boolean showDeleted)</code>
   * method which optionally allows access to deleted paths. Direct access to a deleted contig continues to be valid.
   * It is permissable to delete a contig that is already deleted.
   * This also marks as deleted any path that contains this contig.
   * @param contigId directional identifier of contig being checked.
   * @throws IllegalArgumentException if <code>abs(contigId)</code> is &le; 0 or &gt; <code>numberContigs()</code>.
   */
  void deleteContig(long contigId) throws IllegalArgumentException;

  /**
   * Set the value of an attribute for a contig.
   * @param contigId directional identifier of contig being checked.
   * @param attribute key for an attribute (may be null and should occur in <code>contigAttributes()</code>).
   * @param value the attribute value (may be null).
   * @throws IllegalArgumentException if <code>contigId</code> is &le; 0 or &gt; <code>numberContigs()</code> or has been deleted
   * or attribute is not in <code>contigAttributes()</code>.
   */
  void setContigAttribute(long contigId, String attribute, String value) throws IllegalArgumentException;

  /**
   * Add another path to the graph copying its contents from <code>path</code> and allocating and returning a new unique contig identifier.
   * @param path the path used to set the contigs in the path.
   * @return the identifier for the path (&gt; 0).
   * @throws RuntimeException if <code>path</code> is null.
   */
  long addPath(Path path);

  /**
   * Mark this path as deleted. However, see the <code>paths(long contigId, PathIterator iter, boolean showDeleted)</code>
   * method which optionally allows access to deleted paths. Direct access to a deleted path continues to be valid.
   * It is permissable to delete a path that is already deleted.
   * This has no effect on the contigs in the path.
   * @param pathId directional identifier of path.
   * @throws IllegalArgumentException if <code>abs(pathId)</code> is 0 or &gt; <code>numberPaths()</code>.
   */
  void deletePath(long pathId) throws IllegalArgumentException;

  /**
   * Set the value of an attribute for a path.
   * @param pathId directional identifier of path.
   * @param attribute key for an attribute (may be null and must occur in <code>pathAttributes()</code>).
   * @param value the attribute value (may be null).
   * @throws IllegalArgumentException if <code>abs(pathId)</code> is 0 or &gt; <code>numberPaths()</code> or has been deleted
   * or attribute is not in <code>pathAttributes()</code>.
   */
  void setPathAttribute(long pathId, String attribute, String value);

  /**
   * Create a compacted version of some implementation of the graph with all deleted contigs and paths removed.
   * Contig and path identifiers are not guaranteed to be the same in the new graph.
   * @return a compacted version of this graph.
   */
  MutableGraph compact();

  /**
   * Copy a compacted version of the graph into the supplied parameter
   * @param newGraph destination for the compacted version
   */
  void compact(MutableGraph newGraph);

  /**
   * Add a new attribute to the paths
   * @param name attribute key for an attribute
   * @param description description of the attribute
   */
  void addPathAttribute(String name, String description);

  /**
   * Add a new attribute to the contig
   * @param name attribute key for an attribute
   * @param description description of the attribute
   */
  void addContigAttribute(String name, String description);
}

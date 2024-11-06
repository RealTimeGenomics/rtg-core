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
  void deleteContig(long contigId);

  /**
   * Set the value of an attribute for a contig.
   * @param contigId directional identifier of contig being checked.
   * @param attribute key for an attribute (may be null and should occur in <code>contigAttributes()</code>).
   * @param value the attribute value (may be null).
   * @throws IllegalArgumentException if <code>contigId</code> is &le; 0 or &gt; <code>numberContigs()</code> or has been deleted
   * or attribute is not in <code>contigAttributes()</code>.
   */
  void setContigAttribute(long contigId, String attribute, String value);

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
  void deletePath(long pathId);

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

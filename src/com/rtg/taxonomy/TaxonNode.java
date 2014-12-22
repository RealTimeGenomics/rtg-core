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

package com.rtg.taxonomy;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

/**
 * RTG taxonomy tree node.
 */
public final class TaxonNode implements Comparable<TaxonNode> {

  private final int mId;
  private String mName = null;
  private String mRank = null;

  private TaxonNode mParent = null;

  final TreeSet<TaxonNode> mChildren = new TreeSet<>();

  TaxonNode(int id, String name, String rank) {
    mId = id;
    setName(name);
    setRank(rank);
  }

  TaxonNode(int id) {
    this(id, null, null);
  }

  public int getId() {
    return mId;
  }

  void addChild(TaxonNode node) {
    if (node != null) {
      mChildren.add(node);
      node.mParent = this;
    }
  }

  void detach() {
    if (mParent != null) {
      mParent.mChildren.remove(this);
      mParent = null;
    }
  }

  void setRank(String rank) {
    mRank = rank == null ? null : rank.replaceAll("\\s", " ");
  }

  public String getRank() {
    return mRank;
  }

  void setName(String name) {
    mName = name;
  }

  public String getName() {
    return mName;
  }

  public int getParentId() {
    return mParent == null ? -1 : mParent.mId;
  }

  public TaxonNode getParent() {
    return mParent;
  }

  /**
   * Get a list of the taxon nodes, ordered depth first
   * @return the list of nodes
   */
  public List<TaxonNode> depthFirstTraversal() {
    final ArrayList<TaxonNode> nodes = new ArrayList<>();
    depthFirstTraversal(nodes);
    return nodes;
  }


  /**
   * @return true if this node has no children
   */
  public boolean isLeaf() {
    return mChildren.isEmpty();
  }

  /**
   * @return the number of children
   */
  public int numChildren() {
    return mChildren.size();
  }

  /**
   * @return list of immediate child nodes
   */
  public List<TaxonNode> getChildren() {
    return new ArrayList<>(mChildren);
  }

  private void depthFirstTraversal(ArrayList<TaxonNode> nodes) {
    nodes.add(this);
    for (final TaxonNode child : mChildren) {
      child.depthFirstTraversal(nodes);
    }
  }

  @Override
  public String toString() {
    return getId() + "\t" + getParentId() + "\t" + getRank() + "\t" + getName();
  }

  @Override
  public boolean equals(Object other) {
    return other != null && mId == ((TaxonNode) other).mId;
  }

  @Override
  public int hashCode() {
    return mId;
  }

  @Override
  public int compareTo(TaxonNode o) {
    return mId - o.mId;
  }

}

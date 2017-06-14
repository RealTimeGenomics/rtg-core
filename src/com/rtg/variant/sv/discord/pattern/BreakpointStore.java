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

package com.rtg.variant.sv.discord.pattern;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import com.rtg.util.AutoAddMap;
import com.rtg.util.StringUtils;

/**
*         Date: 16/03/12
*         Time: 10:11 AM
*/
class BreakpointStore implements Iterable<VcfBreakpoint> {

  // Java has no typdef and generic map names are clunky
  static class PositionMap extends TreeSet<VcfBreakpoint> {
  }
  static class RemoteMap extends AutoAddMap<String, PositionMap> {
    @Override
    public PositionMap make() {
      return new PositionMap();
    }
  }
  static class ChrMap extends AutoAddMap<String, RemoteMap> {
    @Override
    public RemoteMap make() {
      return new RemoteMap();
    }
  }

  ChrMap mMap = new ChrMap();

  void add(VcfBreakpoint br) {
    mMap.getOrAdd(br.getLocalChr()).getOrAdd(br.getRemoteChr()).add(br);
  }
  List<String> getChromosomes() {
    final List<String> names = new ArrayList<>();
    names.addAll(mMap.keySet());
    Collections.sort(names);
    return names;
  }
  ChrMap getMap() {
    return mMap;
  }
  @Override
  public Iterator<VcfBreakpoint> iterator() {
    return new StoreIterator(mMap);
  }
  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    forEach(b-> sb.append(b.toString()).append(StringUtils.LS));
    return sb.toString();
  }
  private static class StoreIterator implements Iterator<VcfBreakpoint> {
    Iterator<RemoteMap> mChromosome;
    Iterator<PositionMap> mRemote;
    Iterator<VcfBreakpoint> mPosition;
    StoreIterator(ChrMap map) {
      mChromosome = map.values().iterator();
      assert mChromosome.hasNext();
      mRemote = mChromosome.next().values().iterator();
      assert mRemote.hasNext();
      mPosition = mRemote.next().iterator();
    }

    @Override
    public boolean hasNext() {
      return mChromosome.hasNext() || mRemote.hasNext() || mPosition.hasNext();
    }

    @Override
    public VcfBreakpoint next() {
      while (true) {
        if (!hasNext()) {
          throw new IllegalStateException();
        }
        if (mPosition != null && mPosition.hasNext()) {
          return mPosition.next();
        }
        if (mRemote != null && mRemote.hasNext()) {
          mPosition = mRemote.next().iterator();
          continue;
        }
        if (mChromosome != null && mChromosome.hasNext()) {
          mRemote = mChromosome.next().values().iterator();
        }
      }
    }

    @Override
    public void remove() {
      throw new UnsupportedOperationException();
    }
  }
}

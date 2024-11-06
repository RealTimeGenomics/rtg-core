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

package com.rtg.variant.sv.discord.pattern;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import com.rtg.util.AutoAddMap;
import com.rtg.util.StringUtils;

/**
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
    final List<String> names = new ArrayList<>(mMap.keySet());
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
    forEach(b-> sb.append(b).append(StringUtils.LS));
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

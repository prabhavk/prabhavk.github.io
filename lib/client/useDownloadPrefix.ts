"use client";

import { useEffect, useState } from "react";

const KEY = "emtr:downloadPrefix";

export function useDownloadPrefix() {
  const [prefix, setPrefix] = useState<string>("");

  useEffect(() => {
    try {
      const saved = localStorage.getItem(KEY);
      if (saved !== null) setPrefix(saved);
    } catch {}
  }, []);

  useEffect(() => {
    try { localStorage.setItem(KEY, prefix); } catch {}
  }, [prefix]);

  return [prefix, setPrefix] as const;
}

export function getDownloadPrefix(): string {
  if (typeof window === "undefined") return "";
  return localStorage.getItem(KEY) ?? "";
}

/** Utility to make safe filenames. */
export function makeDownloadName(fallbackBase: string, ext?: string): string {
  const base = (getDownloadPrefix().trim() || fallbackBase).replace(/\s+/g, "-");
  const safe = base.replace(/[^\w.\-]+/g, "_");
  return ext ? `${safe}.${ext}` : safe;
}

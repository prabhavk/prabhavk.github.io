// components/rep-context.tsx
"use client";
import React, { createContext, useContext, useEffect, useState } from "react";

type RepCtx = { rep: string | null; setRep: (r: string | null) => void };
const Ctx = createContext<RepCtx>({ rep: null, setRep: () => {} });

const keyFor = (job: string) => `emtr:selectedRep:${job}`;

export function RepetitionProvider({ job, children }: { job: string; children: React.ReactNode }) {
  const [rep, setRep] = useState<string | null>(null);
  useEffect(() => {
    if (!job) return;
    try { setRep(localStorage.getItem(keyFor(job))); } catch {}
  }, [job]);
  useEffect(() => {
    if (!job) return;
    try {
      rep ? localStorage.setItem(keyFor(job), rep) : localStorage.removeItem(keyFor(job));
    } catch {}
  }, [job, rep]);
  return <Ctx.Provider value={{ rep, setRep }}>{children}</Ctx.Provider>;
}

export const useRep = () => useContext(Ctx);

"use client";

import React, {
  createContext,
  useContext,
  useLayoutEffect,
  useRef,
  useState,
  useCallback,
} from "react";

type Status = "idle" | "started" | "done" | "error";

type ProgressCtx = {
  jobId?: string;
  status: Status;
  logs: string[];
  start: (jobId: string) => void;
  append: (line: string) => void;
  stop: () => void;
  clear: () => void;
  attachScroller: (el: HTMLElement | null) => void;
};

const Ctx = createContext<ProgressCtx | null>(null);

export function useProgress() {
  const ctx = useContext(Ctx);
  if (!ctx) throw new Error("useProgress must be used within <ProgressProvider>");
  return ctx;
}

export default function ProgressProvider({ children }: { children: React.ReactNode }) {
  const [jobId, setJobId] = useState<string | undefined>(undefined);
  const [status, setStatus] = useState<Status>("idle");
  const [logs, setLogs] = useState<string[]>([]);
  const pollRef = useRef<number | null>(null);

  // element to auto-scroll
  const scrollerRef = useRef<HTMLElement | null>(null);
  const attachScroller = useCallback((el: HTMLElement | null) => {
    scrollerRef.current = el ?? null;
  }, []);

  const clear = useCallback(() => {
    setLogs([]);
    setStatus("idle");
    setJobId(undefined);
    if (pollRef.current) window.clearInterval(pollRef.current);
    pollRef.current = null;
  }, []);

  const append = useCallback((line: string) => {
    setLogs((prev) => [...prev, line]);
  }, []);

  // Auto-scroll to bottom whenever logs/status change
  useLayoutEffect(() => {
    const el = scrollerRef.current;
    if (!el) return;

    const nearBottom = el.scrollHeight - el.scrollTop - el.clientHeight < 40;

    if (nearBottom || status === "started") {
      // wait for layout/paint to account for the new line
      requestAnimationFrame(() => {
        const prev = el.style.scrollBehavior;
        // smooth scrolling only while actively "started"
        el.style.scrollBehavior = status === "started" ? "smooth" : "auto";
        el.scrollTop = el.scrollHeight;
        el.style.scrollBehavior = prev || "";
      });
    }
  }, [logs, status]);

  // Call this when you kick off EMTR
  const start = useCallback((id: string) => {
    setJobId(id);
    setStatus("started");
    setLogs([]);
    if (pollRef.current) {
      window.clearInterval(pollRef.current);
      pollRef.current = null;
    }
  }, []);

  // Call this when the worker finishes or errors
  const stop = useCallback(() => {
    setStatus("idle"); // or set to "done"/"error" if you want transient states
    if (pollRef.current) {
      window.clearInterval(pollRef.current);
      pollRef.current = null;
    }
  }, []);

  return (
    <Ctx.Provider
      value={{ jobId, status, logs, start, append, stop, clear, attachScroller }}
    >
      {children}
    </Ctx.Provider>
  );
}

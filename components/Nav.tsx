// components/Nav.tsx
'use client';

import Link from "next/link";
import { usePathname } from "next/navigation";

type Tab = { href: string; label: string; mode: "exact" | "deep" };

const tabs: Tab[] = [
  { href: "/",        label: "Input Form",             mode: "exact" },
  { href: "/runs",    label: "Precomputed Results",    mode: "deep"  },
  { href: "/violin",  label: "MLL D,P&S",              mode: "deep"  },
  { href: "/mle",     label: "MLE in D",               mode: "deep"  },
  { href: "/emconv",  label: "ECDLL",                    mode: "deep"  },
  { href: "/wmw",     label: "Wilcoxon-Mann-Whitney",  mode: "deep"  },  
  { href: "/sigrec",  label: "Significance of Recall", mode: "deep"  },  
  { href: "/spiral",  label: "Spiral Plots",           mode: "deep"  },  
  
];

export default function Nav() {
  const rawPath = usePathname() ?? "/";
  const pathname = rawPath !== "/" && rawPath.endsWith("/") ? rawPath.slice(0, -1) : rawPath;

  const linkBase =
    "inline-flex items-center justify-center text-center px-4 py-2 rounded-md text-sm font-semibold " +
    "transition-colors outline-none focus-visible:ring-2 focus-visible:ring-yellow-400 " +
    "focus-visible:ring-offset-2 focus-visible:ring-offset-black";

  // const active = "bg-gray-200 text-black shadow-md";
  const active    = "bg-gray-200 text-black shadow-md hover:bg-white hover:!text-black";
  const idleDefault = "bg-gray-700 text-white hover:bg-yellow-500";
  const idleRuns    = "bg-gray-700 text-white hover:bg-white hover:!text-black";

  const isActive = (t: Tab) =>
    t.mode === "exact" ? pathname === t.href : pathname === t.href || pathname.startsWith(t.href + "/");

  return (
    <nav className="w-full border-b bg-black">
      <div className="max-w-6xl mx-auto px-4 py-3 flex flex-wrap gap-3 justify-center">
        {tabs.map((t) => {
          const on = isActive(t);
          const idle = t.href === "/runs" ? idleRuns : idleDefault;
          return (
            <Link
              key={t.href}
              href={t.href}
              className={`${linkBase} ${on ? active : idle}`}
              aria-current={on ? "page" : undefined}
            >
              {t.label}
            </Link>
          );
        })}
      </div>
    </nav>
  );
}

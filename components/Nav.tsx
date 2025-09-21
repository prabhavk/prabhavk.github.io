// components/Nav.tsx
'use client';

import Link from "next/link";
import { usePathname } from "next/navigation";

type Tab = { href: string; label: string; mode: "exact" | "deep" };

const tabs: Tab[] = [
  { href: "/",                label: "Input",                 mode: "exact" }, // replace progress pane with three-stained clock
  { href: "/runs",            label: "Available Results",     mode: "deep"  },  
  { href: "/cake",            label: "Cake",                  mode: "deep"  },
  { href: "/violin_methods",  label: "Violin methods",        mode: "deep"  }, // pack in same page and add layers
  { href: "/violin_roots",    label: "Violin roots",          mode: "deep"  }, // pack in same page and add layers
  { href: "/ecdll",           label: "ECDLL",                 mode: "deep"  }, // add layers 
  { href: "/wmw_method",      label: "WMW method",            mode: "deep"  }, // pack in same page and add layers
  { href: "/wmw_root",        label: "WMW root",              mode: "deep"  }, // pack in same page and add layers
  { href: "/mle    ",         label: "MLE",                   mode: "deep"  }, // add layers 
];

export default function Nav() {
  const rawPath = usePathname() ?? "/";
  const pathname = rawPath !== "/" && rawPath.endsWith("/") ? rawPath.slice(0, -1) : rawPath;

  const isActive = (t: Tab) =>
    t.mode === "exact" ? pathname === t.href : pathname === t.href || pathname.startsWith(t.href + "/");

  const linkBase =
    "inline-flex items-center justify-center text-center px-4 py-2 rounded-md text-sm font-semibold " +
    "transition-colors outline-none focus-visible:ring-2 focus-visible:ring-green-500 " +
    "focus-visible:ring-offset-2 focus-visible:ring-offset-black text-black";

  // selected red; idle green; both hover yellow; text stays black
  const active      = "bg-gray-300 hover:bg-green-600";
  const idleGeneral = "bg-gray-500 hover:bg-blue-600";

 
  return (
    <nav className="w-full border-b bg-[#3C1421]">
      {/* Top row: brand only */}
      <div className="max-w-6xl mx-auto px-4 py-2 flex items-center">        
      </div>

      {/* Tabs row */}
      <div className="max-w-6xl mx-auto px-4 pb-3 flex flex-wrap gap-3 justify-center">
        {tabs.map((t) => {
          const on = isActive(t);
          return (
            <Link
              key={t.href}
              href={t.href}
              className={`${linkBase} ${on ? active : idleGeneral}`}
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

"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";

export default function Navbar() {
  const pathname = usePathname();

  const linkClasses = (path: string) =>
    `px-3 py-2 rounded ${
      pathname === path
        ? "bg-gray-600 font-semibold text-white"
        : "hover:bg-black text-white"
    }`;

  return (
    <nav className="bg-red-600 px-6 py-3 flex gap-4">
      <Link href="/" className={linkClasses("/")}>Form</Link>
      <Link href="/results" className={linkClasses("/results")}>Results</Link>
    </nav>
  );
}

import "./globals.css";
import Link from "next/link";

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en">
      <body className="bg-gray-50 text-black">
        {/* Navbar */}
        <nav className="bg-blue-600 text-white px-6 py-3 flex space-x-6">
          <Link href="/" className="hover:underline">
            Form
          </Link>
          <Link href="/results" className="hover:underline">
            Results
          </Link>
        </nav>
        <div className="p-6">{children}</div>
      </body>
    </html>
  );
}

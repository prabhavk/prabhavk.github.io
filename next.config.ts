// next.config.ts
const nextConfig = {
  output: "export",
  ...(process.env.NODE_ENV === "development"
    ? {
        experimental: {
          allowedDevOrigins: ["http://192.168.29.107:3000"],
        },
      }
    : {}),
};

export default nextConfig;

// next.config.ts
const nextConfig = {
  // output: "export",
  allowedDevOrigins: ["http://192.168.29.107:3000"],
  ...(process.env.NODE_ENV === "development"
    ? {
        experimental: {          
        },
      }
    : {}),
};

export default nextConfig;

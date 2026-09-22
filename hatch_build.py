import os

from hatchling.metadata.plugin.interface import MetadataHookInterface


class CustomMetadataHook(MetadataHookInterface):
    """Override wheel metadata for GCP bleeding-edge development builds.

    This overrides the version according to:
    <version>.dev<BITBUCKET_BUILD_NUMBER>+g<BITBUCKET_COMMIT:0:8>
    """

    def update(self, metadata: dict) -> None:
        build_version = os.environ.get("PYLINAC_BUILD_VERSION")
        if build_version is not None:
            metadata["version"] = build_version
